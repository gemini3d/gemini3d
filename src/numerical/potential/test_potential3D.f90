program test_potential3D

use mpi

implicit none (type, external)

type(mumps_struc) :: mumps_par

integer :: ierr

integer, parameter :: npts1=256,npts2=256,npts3=12
integer, parameter :: lk=npts1*npts2*npts3
integer :: lent
integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(8), dimension(:), allocatable :: M
real(8), dimension(:), allocatable :: b
real(8) :: dx1
real(8), dimension(npts2,npts3) :: Vminx1,Vmaxx1
real(8), dimension(npts1,npts3) :: Vminx2,Vmaxx2
real(8), dimension(npts1,npts2) :: Vminx3,Vmaxx3
real(8), dimension(:,:), allocatable ::  Mfull
real(8) :: tstart,tfin


!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
lent=7*(npts1-2)*(npts2-2)*(npts3-2)                                           !interior entries
lent=lent+2*(npts1-2)*(npts2-2)+2*(npts2-2)*(npts3-2)+2*(npts1-2)*(npts3-2)    !6 faces of cube
lent=lent+4*(npts1-2)+4*(npts2-2)+4*(npts3-2)                                  !12 edges
lent=lent+8                                                                    !8 corners
allocate(ir(lent),ic(lent),M(lent),b(lk))
lx1=npts1
lx2=npts2
lx3=npts3

dx1=1.0/npts1           !scale dx so the domain of problem is [0,1]

Vminx1(:,:)=0
Vmaxx1(:,:)=0
Vminx2(:,:)=0
Vmaxx2(:,:)=0
Vminx3(:,:)=0
Vmaxx3(:,:)=10

M(:)=0.0
b(:)=0.0
ient=1


!LOAD UP MATRIX ELEMENTS
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1     !linear index referencing Phi(ix1,ix2,ix3) as a column vector.  Also row # of big matrix

      if (ix1==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx1(ix2,ix3)
        ient=ient+1
      elseif (ix1==lx1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx1(ix2,ix3)
        ient=ient+1
      elseif (ix2==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx2(ix1,ix3)
        ient=ient+1
      elseif (ix2==lx2) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx2(ix1,ix3)
        ient=ient+1
      elseif (ix3==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx3(ix1,ix2)
        ient=ient+1
      elseif (ix3==lx3) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx3(ix1,ix2)
        ient=ient+1
      else                       !INTERIOR
        !ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1*lx2
        M(ient)=1.0
        ient=ient+1

        !ix1,ix2-1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=1.0
        ient=ient+1

        !ix1-1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=1.0
        ient=ient+1

        !ix1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-6.0
        ient=ient+1

        !ix1+1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=1.0
        ient=ient+1

        !ix1,ix2+1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=1.0
        ient=ient+1

        !ix1,ix2,ix3+1
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1*lx2
        M(ient)=1.0
        ient=ient+1
      end if
    end do
  end do
end do


!CORRECT FOR DX /= 1
b=b*dx1**2


!OUTPUT FULL MATRIX FOR DEBUGGING IF ITS NOT TOO BIG (ZZZ --> CAN BE COMMENTED OUT)
block
  integer :: u
  open(newunit=u,file='test_potential3D.dat',status='replace')
  write(u,*) lx1,lx2,lx3
  if (lk<150) then
    allocate(Mfull(lk,lk))
    Mfull(:,:)=0.0
    do ient=1,size(ir)
      Mfull(ir(ient),ic(ient))=M(ient)
    end do
    call write2Darray(u,Mfull)
    call writearray(u,b)
    deallocate(Mfull)
  end if


  !------------------------------------------------------------
  !-------DO SOME STUFF TO CALL MUMPS
  !------------------------------------------------------------
  call MPI_INIT(IERR)
  if (ierr/=0) error stop 'mpi init'

  ! Define a communicator for the package.
  mumps_par%COMM = MPI_COMM_WORLD


  !Initialize an instance of the package
  !for L U factorization (sym = 0, with working host)
  mumps_par%JOB = -1
  mumps_par%SYM = 0
  mumps_par%PAR = 1
  call DMUMPS(mumps_par)


  !Define problem on the host (processor 0)
  if ( mumps_par%MYID .eq. 0 ) then
    mumps_par%N=lk
    mumps_par%NZ=lent
    allocate( mumps_par%IRN ( mumps_par%NZ ) )
    allocate( mumps_par%JCN ( mumps_par%NZ ) )
    allocate( mumps_par%A( mumps_par%NZ ) )
    allocate( mumps_par%RHS ( mumps_par%N  ) )
    mumps_par%IRN=ir
    mumps_par%JCN=ic
    mumps_par%A=M
    mumps_par%RHS=b

  !  mumps_par%ICNTL(7)=6    !force a particular reordering - see mumps docs
  !  mumps_par%ICNTL(28)=2
  !  mumps_par%ICNTL(29)=2
  end if


  !Call package for solution
  mumps_par%JOB = 6
  call cpu_time(tstart)
  call DMUMPS(mumps_par)
  call cpu_time(tfin)
  write(*,*) 'Solve took ',tfin-tstart,' seconds...'


  !Solution has been assembled on the host
  if ( mumps_par%MYID == 0 ) then
    call writearray(u,mumps_par%RHS/dx1**2)
  end if
  close(u)
end block

!Deallocate user data
if ( mumps_par%MYID == 0 ) then
  deallocate( mumps_par%IRN )
  deallocate( mumps_par%JCN )
  deallocate( mumps_par%A   )
  deallocate( mumps_par%RHS )
end if
deallocate(ir,ic,M,b)


!Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
call DMUMPS(mumps_par)

call MPI_FINALIZE(IERR)
if (ierr /= 0) error stop 'mpi finalize



contains

  subroutine writearray(fileunit,array)
    integer, intent(in) :: fileunit
    real(8), dimension(:), intent(in) :: array

    integer :: k

    do k=1,size(array)
      write(fileunit,*) array(k)
    end do
  end subroutine writearray


  subroutine write2Darray(fileunit,array)
    integer, intent(in) :: fileunit
    real(8), dimension(:,:), intent(in) :: array

    integer :: k1,k2

    do k1=1,size(array,1)
      write(fileunit,'(f4.0)') (array(k1,k2), k2=1,size(array,2))
    end do
  end subroutine write2Darray

end program test_potential3D
