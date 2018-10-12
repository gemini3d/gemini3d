module potential

!use calculus
implicit none
include 'mpif.h'
include 'dmumps_struc.h'

integer, dimension(:), pointer, save :: mumps_perm   !cached permutation, unclear whether save is necessary...

contains


  function elliptic2D(srcterm,SigP,SigH,Vminx2,Vmaxx2,Vminx3,Vmaxx3,dx2,dx2i,dx3,dx3i,perflag)

    !------------------------------------------------------------
    !-------SOLVE POISSONS'S EQUATION IN 2D USING MUMPS
    !------------------------------------------------------------

    real(wp), dimension(:,:), intent(in) :: srcterm,SigP,SigH
    real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
    real(wp), dimension(:), intent(in) :: Vminx3,Vmaxx3
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(:), intent(in) :: dx2i
    real(wp), dimension(0:), intent(in) :: dx3
    real(wp), dimension(:), intent(in) :: dx3i
    logical, intent(in) :: perflag

    real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh2
    real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh3
    real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: gradSigH2,gradSigH3

    integer :: ix2,ix3,lx2,lx3
    integer :: lPhi,lent
    integer :: iPhi,ient
    integer, dimension(:), allocatable :: ir,ic
    real(wp), dimension(:), allocatable :: M
    real(wp), dimension(:), allocatable :: b
    real(wp) :: tstart,tfin
    type (DMUMPS_STRUC) mumps_par
    integer :: myid, ierr

    real(wp), dimension(size(SigP,1),size(SigP,2)) :: elliptic2D


    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)


    !ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
    if (myid==0) then
      lx2=size(SigP,1)
      lx3=size(SigP,2)
      lPhi=lx2*lx3
      lent=5*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)
      allocate(ir(lent),ic(lent),M(lent),b(lPhi))


      !PREP INPUT DATA FOR SOLUTION OF SYSTEM
      SigPh2(1,:)=0.0
      SigPh2(2:lx2,:)=0.5*(SigP(1:lx2-1,:)+SigP(2:lx2,:))
      SigPh3(:,1)=0.0
      SigPh3(:,2:lx3)=0.5*(SigP(:,1:lx3-1)+SigP(:,2:lx3))

!      gradSigH2=grad2D1(SigH,dx2(1:lx2))
!      gradSigH3=grad2D1(SigH,dx3(1:lx3))
      gradSigH2=0d0
      gradSigH3=0d0

      !------------------------------------------------------------
      !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
      !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
      !-------GUIDE).
      !------------------------------------------------------------

      !LOAD UP MATRIX ELEMENTS
      M(:)=0.0
      b=pack(srcterm,.true.)           !boundaries overwritten later
      ient=1
      do ix3=1,lx3
        do ix2=1,lx2
          iPhi=lx2*(ix3-1)+ix2     !linear index referencing Phi(ix2,ix3) as a column vector.  Also row of big matrix

          if (ix2==1) then          !BOTTOM GRID POINTS + CORNER
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1.0
            b(iPhi)=Vminx2(ix3)
            ient=ient+1
          elseif (ix2==lx2) then    !TOP GRID POINTS + CORNER
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1.0
            b(iPhi)=Vmaxx2(ix3)
            ient=ient+1
          elseif (ix3==1) then      !LEFT BOUNDARY
            ir(ient)=iPhi
            ic(ient)=iPhi  
            M(ient)=1.0
            b(iPhi)=Vminx3(ix2)
            ient=ient+1
          elseif (ix3==lx3) then    !RIGHT BOUNDARY
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=1.0
            b(iPhi)=Vmaxx3(ix2)
            ient=ient+1
          else                      !INTERIOR
            !ix2,ix3-1 grid point in ix2,ix3 equation
            ir(ient)=iPhi
            ic(ient)=iPhi-lx2
            M(ient)=SigPh3(ix2,ix3)/(dx3i(ix3)*dx3(ix3))+gradSigH2(ix2,ix3)/(dx3(ix3)+dx3(ix3+1))
            ient=ient+1

            !ix2-1,ix3 grid point
            ir(ient)=iPhi
            ic(ient)=iPhi-1
            M(ient)=SigPh2(ix2,ix3)/(dx2i(ix2)*dx2(ix2))-gradSigH3(ix2,ix3)/(dx2(ix2)+dx2(ix2+1))
            ient=ient+1

            !ix2,ix3 grid point
            ir(ient)=iPhi
            ic(ient)=iPhi
            M(ient)=-1.0*SigPh2(ix2+1,ix3)/(dx2i(ix2)*dx2(ix2+1)) &
                    -1.0*SigPh2(ix2,ix3)/(dx2i(ix2)*dx2(ix2)) &
                    -1.0*SigPh3(ix2,ix3+1)/(dx3i(ix3)*dx3(ix3+1)) &
                    -1.0*SigPh3(ix2,ix3)/(dx3i(ix3)*dx3(ix3))
            ient=ient+1

            !ix2+1,ix3 grid point
            ir(ient)=iPhi
            ic(ient)=iPhi+1
            M(ient)=SigPh2(ix2+1,ix3)/(dx2i(ix2)*dx2(ix2+1))+gradSigH3(ix2,ix3)/(dx2(ix2)+dx2(ix2+1))
            ient=ient+1

            !ix2,ix3+1 grid point
            ir(ient)=iPhi
            ic(ient)=iPhi+lx2
            M(ient)=SigPh3(ix2,ix3+1)/(dx3i(ix3)*dx3(ix3+1))-gradSigH2(ix2,ix3)/(dx3(ix3)+dx3(ix3+1))
            ient=ient+1
          end if
        end do
      end do
    end if


    !FIRE UP MUMPS
    mumps_par%COMM = MPI_COMM_WORLD
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 1
    call DMUMPS(mumps_par)


    !LOAD OUR PROBLEM
    if ( mumps_par%MYID==0 ) then
      mumps_par%N=lPhi
      mumps_par%NZ=lent
      allocate( mumps_par%IRN ( mumps_par%NZ ) )
      allocate( mumps_par%JCN ( mumps_par%NZ ) )
      allocate( mumps_par%A( mumps_par%NZ ) )
      allocate( mumps_par%RHS ( mumps_par%N  ) )
      mumps_par%IRN=ir
      mumps_par%JCN=ic
      mumps_par%A=M
      mumps_par%RHS=b
      deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

      if (perflag) then       !used cached permutation
        allocate(mumps_par%PERM_IN(mumps_par%N))
        mumps_par%PERM_IN=mumps_perm
        mumps_par%ICNTL(7)=1
      end if
    end if


    !SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
    mumps_par%JOB = 6
    call DMUMPS(mumps_par)


    !STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
    !(can save ~25% execution time and improves scaling with openmpi
    ! ~25% more going from 1-2 processors)
    if ( mumps_par%MYID==0 ) then
      mumps_perm=mumps_par%SYM_PERM
      elliptic2D=reshape(mumps_par%RHS,[lx2,lx3])

      deallocate( mumps_par%IRN )
      deallocate( mumps_par%JCN )
      deallocate( mumps_par%A   )
      deallocate( mumps_par%RHS )
    end if
    mumps_par%JOB = -2
    call DMUMPS(mumps_par)
  end function elliptic2D


  function poisson2D(rho,Vminx1,Vmaxx1,Vminx2,Vmaxx2,dx1,perflag)

    !------------------------------------------------------------
    !-------SOLVE POISSONS'S EQUATION IN 2D USING MUMPS
    !------------------------------------------------------------

    real(wp), dimension(:,:), intent(in) :: rho
    real(wp), dimension(:), intent(in) :: Vminx1,Vmaxx1
    real(wp), dimension(:), intent(in) :: Vminx2,Vmaxx2
    real(wp), intent(in) :: dx1
!    real(wp), dimension(:), intent(in) :: dx2
    logical, intent(in) :: perflag

    integer :: ix1,ix2,lx1,lx2
    integer :: lPhi, lent
    integer :: iPhi,ient
    integer, dimension(:), allocatable :: ir,ic
    real(wp), dimension(:), allocatable :: M
    real(wp), dimension(:), allocatable :: b
    real(wp) :: tstart,tfin
    type (DMUMPS_STRUC) mumps_par
    integer :: myid, ierr
    real(wp), dimension(size(rho,1),size(rho,2)) :: poisson2D


    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

    if (myid==0) then
    !------------------------------------------------------------
    !-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
    !-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
    !-------GUIDE).
    !------------------------------------------------------------
    lx1=size(rho,1)
    lx2=size(rho,2)
    lPhi=lx1*lx2
    lent=5*(lx1-2)*(lx2-2)+2*lx1+2*(lx2-2)
    allocate(ir(lent),ic(lent),M(lent),b(lPhi))


    !LOAD UP MATRIX ELEMENTS
    M(:)=0.0
    b=pack(rho,.true.)           !boundaries overwritten later
    ient=1
    do ix2=1,lx2
      do ix1=1,lx1
        iPhi=lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2) as a column vector.  Also row of big matrix

        if (ix1==1) then          !BOTTOM GRID POINTS + CORNER
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vminx1(ix2)
          ient=ient+1
        elseif (ix1==lx1) then    !TOP GRID POINTS + CORNER
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vmaxx1(ix2)
          ient=ient+1
        elseif (ix2==1) then      !LEFT BOUNDARY
          ir(ient)=iPhi
          ic(ient)=iPhi  
          M(ient)=1.0
          b(iPhi)=Vminx2(ix1)
          ient=ient+1
        elseif (ix2==lx2) then    !RIGHT BOUNDARY
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1.0
          b(iPhi)=Vmaxx2(ix1)
          ient=ient+1
        else                      !INTERIOR
          !ix1,ix2-1 grid point in ix1,ix2 equation
          ir(ient)=iPhi
          ic(ient)=iPhi-lx1
          M(ient)=1.0
          ient=ient+1

          !ix1-1,ix2 grid point
          ir(ient)=iPhi
          ic(ient)=iPhi-1
          M(ient)=1.0
          ient=ient+1

          !ix1,ix2 grid point
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=-4.0
          ient=ient+1

          !ix1+1,ix2 grid point
          ir(ient)=iPhi
          ic(ient)=iPhi+1
          M(ient)=1.0
          ient=ient+1

          !ix1,ix2+1 grid point
          ir(ient)=iPhi
          ic(ient)=iPhi+lx1
          M(ient)=1.0
          ient=ient+1
        end if
      end do
    end do


    !CORRECT FOR DX /= 1
    b=b*dx1**2
    end if

    !FIRE UP MUMPS
    mumps_par%COMM = MPI_COMM_WORLD
    mumps_par%JOB = -1
    mumps_par%SYM = 0
    mumps_par%PAR = 1
    call DMUMPS(mumps_par)


    !LOAD OUR PROBLEM
    if ( mumps_par%MYID .eq. 0 ) then
      mumps_par%N=lPhi
      mumps_par%NZ=lent
      allocate( mumps_par%IRN ( mumps_par%NZ ) )
      allocate( mumps_par%JCN ( mumps_par%NZ ) )
      allocate( mumps_par%A( mumps_par%NZ ) )
      allocate( mumps_par%RHS ( mumps_par%N  ) )
      mumps_par%IRN=ir
      mumps_par%JCN=ic
      mumps_par%A=M
      mumps_par%RHS=b
      deallocate(ir,ic,M,b)     !clear memory before solve begins!!!

      if (perflag) then       !used cached permutation
        allocate(mumps_par%PERM_IN(mumps_par%N))
        mumps_par%PERM_IN=mumps_perm
        mumps_par%ICNTL(7)=1
      end if
    end if


    !SOLVE OUR PROBLEM (ALL WORKERS NEED TO SEE THIS CALL)
    mumps_par%JOB = 6
    call cpu_time(tstart)
    call DMUMPS(mumps_par)
    call cpu_time(tfin)
    write(*,*) 'Solve took ',tfin-tstart,' seconds...'


    !STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
    !(can save ~25% execution time and improves scaling with openmpi
    ! ~25% more going from 1-2 processors)
    if ( mumps_par%MYID .eq. 0 ) then
      mumps_perm=mumps_par%SYM_PERM
      poisson2D=reshape(mumps_par%RHS/dx1**2,[lx1,lx2])

      deallocate( mumps_par%IRN )
      deallocate( mumps_par%JCN )
      deallocate( mumps_par%A   )
      deallocate( mumps_par%RHS )
    end if
    mumps_par%JOB = -2
    call DMUMPS(mumps_par)
  end function poisson2D

end module potential
