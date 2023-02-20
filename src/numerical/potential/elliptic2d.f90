submodule (PDEelliptic) elliptic2d

implicit none (type, external)

contains

module procedure elliptic2D_polarization

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!-------INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!-------TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!-------THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!-------ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!-------BY ZEROING OUT THE INERTIAL CAPACITANCE.
!-------
!-------THIS FORM IS INTENDED TO  WORK WITH CURVILINEAR MESHES.
!-------NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!-------BE USED HERE!!!
!-------The equation solved by this subroutine is:
!-------
!-------  d/dx2(A dV/dx2) + d/dx3(A' dV/dx3) + B dV/dx2 - C dV/x3 + ...
!-------  d/dx2(D d/dt(dV/dx2)) + d/dx3(D d/dt(dV/dx3)) + ...
!-------  d/dx2( D*v2 d^V/dx2^2 + D*v3*d^2V/dx3/dx2 ) + ...
!-------  d/dx3( D*v2 d^2V/dx2/dx3 + D*v3 d^2V/dx3^2 ) = srcterm
!-------
!-------  for GEMINI:  A=SigP2, A'=SigP3, B=d/dx3(SigH), C=d/dx2(SigH),
!-------               D=Cm
!------------------------------------------------------------

real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigPh2    !I'm too lazy to recode these as SigP2h2, etc.
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: SigPh3
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: Cmh2
real(wp), dimension(1:size(SigP2,1),1:size(SigP2,2)) :: Cmh3

real(wp) :: coeff    !coefficient for calculating polarization terms
integer :: ix2,ix3,lx2,lx3    !this overwrites the
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b

type(MUMPS_STRUC) :: mumps_par

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then
lx2=size(SigP2,1)    !note that these are full-grid sizes since grid module globals are not in scope
lx3=size(SigP2,2)
lPhi=lx2*lx3
!      lent=5*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)    !static model; left here as a reference
lent=17*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)-3*2*(lx2-2)-3*2*(lx3-2)
!! interior+boundary-x3_adj-x2_adj.  Note that are 3 sets of entries for each adjacent point
allocate(ir(lent),ic(lent),M(lent),b(lPhi))
if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx2,lx3
if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


!PREP INPUT DATA FOR SOLUTION OF SYSTEM
SigPh2(1,:)= 0
SigPh2(2:lx2,:)=0.5_wp*(SigP2(1:lx2-1,:)+SigP2(2:lx2,:))
!! note the different conductiances here to be associated with derivatives in different directions
SigPh3(:,1)= 0
SigPh3(:,2:lx3)=0.5_wp*(SigP3(:,1:lx3-1)+SigP3(:,2:lx3))
Cmh2(1,:)= 0
Cmh2(2:lx2,:)=0.5_wp*(Cm(1:lx2-1,:)+Cm(2:lx2,:))
Cmh3(:,1)= 0
Cmh3(:,2:lx3)=0.5_wp*(Cm(:,1:lx3-1)+Cm(:,2:lx3))


!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
!LOAD UP MATRIX ELEMENTS
M(:)= 0
b = pack(srcterm,.true.)           !boundaries overwritten later, polarization terms also added later.
ient=1

loopx3: do ix3=1,lx3
  loopx2: do ix2=1,lx2
    iPhi=lx2*(ix3-1)+ix2     !linear index referencing Phi(ix2,ix3) as a column vector.  Also row of big matrix

    if (ix2==1) then
      !! BOTTOM GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vminx2(ix3)
      ient=ient+1
      cycle
    elseif (ix2==lx2) then
      !! TOP GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vmaxx2(ix3)
      ient=ient+1
      cycle
    elseif (ix3==1) then
      !! LEFT BOUNDARY
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vminx3(ix2)
      ient=ient+1
      cycle
    elseif (ix3==lx3) then
      !! RIGHT BOUNDARY
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vmaxx3(ix2)
      ient=ient+1
      cycle
    endif

    !! INTERIOR LOCATION
    !> ix2-1,ix3-2 grid point
    !>>  because we are one interior point in for x3, ix3-2 is "before" the boundary; we'll assume it's at the
    !>>  same potential.  Effectively the electric field is assumed to go to zero at the boundary.
    coeff=-Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    if (ix3==2) then    !out of bounds, use nearest BC, and add to known vector
      b(iPhi)=b(iPhi)-coeff*Vminx3(ix2-1)
    else    !in bounds, add to matrix
      ir(ient)=iPhi
      ic(ient)=iPhi-2*lx2-1

      M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

      ient=ient+1
    end if


    !> ix2,ix3-2 grid point
    coeff=-Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)*dx3iall(ix3-1)) )
    if (ix3==2) then
    !! bit of intentional code duplication here and in the following sections to keep things organized in a way I can debug...
      b(iPhi)=b(iPhi)-coeff*Vminx3(ix2)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-2*lx2

      M(ient)=coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

      ient=ient+1
    end if


    !> ix2+1,ix3-2 grid point
    coeff=Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    if (ix3==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx3(ix2+1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-2*lx2+1

      M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

      ient=ient+1
    end if


    !> ix2-2,ix3-1
    coeff=-Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(ix3-1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-lx2-2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !!!ix2,ix3-1 grid point in ix2,ix3 equation
    ir(ient)=iPhi
    ic(ient)=iPhi-lx2

    M(ient)=SigPh3(ix2,ix3)/(dx3iall(ix3)*dx3all(ix3))+gradSigH2(ix2,ix3)/(dx3all(ix3)+dx3all(ix3+1))    !static terms

    coeff=Cmh3(ix2,ix3)/(dt*dx3iall(ix3)*dx3all(ix3))
    M(ient)=M(ient)+coeff   !polarization time derivative terms
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3-1)
    !! add in polarziation terms that include previous time step potential at this grid point

    coeff=Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3)*dx3iall(ix3-1)) )+ &
    Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)*dx3iall(ix3-1)) )
    M(ient)=M(ient)+coeff   !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) )+ &
    Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    M(ient)=M(ient)+coeff   !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

    ient=ient+1


    !> ix2+2,ix3-1 grid point
    coeff=-Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2+1)+dx2all(ix2+2))*(dx2all(ix2)+dx2all(ix2+1))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3-1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-lx2+2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !> ix2-2,ix3 grid point
    coeff=-Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)*dx2iall(ix2-1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(ix3)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-2

      M(ient)=coeff
      !! d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

      ient=ient+1
    end if


    !> ix2-1,ix3 grid point
    ir(ient)=iPhi
    ic(ient)=iPhi-1

    M(ient)=SigPh2(ix2,ix3)/(dx2iall(ix2)*dx2all(ix2))-gradSigH3(ix2,ix3)/(dx2all(ix2)+dx2all(ix2+1))    !static

    coeff=Cmh2(ix2,ix3)/(dt*dx2iall(ix2)*dx2all(ix2))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2-1,ix3)    !BC's and pol. time deriv.

    coeff=Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2)*dx2iall(ix2-1)) )+ &
    Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)*dx2iall(ix2-1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) )+ &
    Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    M(ient)=M(ient)+coeff     !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

    ient=ient+1


    !!!ix2,ix3 grid point (main diagonal)
    ir(ient)=iPhi
    ic(ient)=iPhi

    M(ient)=-SigPh2(ix2+1,ix3)/(dx2iall(ix2)*dx2all(ix2+1)) &
    -SigPh2(ix2,ix3)/(dx2iall(ix2)*dx2all(ix2)) &
    -SigPh3(ix2,ix3+1)/(dx3iall(ix3)*dx3all(ix3+1)) &
    -SigPh3(ix2,ix3)/(dx3iall(ix3)*dx3all(ix3))    !static

    coeff=-Cmh2(ix2+1,ix3)/(dt*dx2iall(ix2)*dx2all(ix2+1)) &
    -Cmh2(ix2,ix3)/(dt*dx2iall(ix2)*dx2all(ix2)) &
    -Cmh3(ix2,ix3+1)/(dt*dx3iall(ix3)*dx3all(ix3+1)) &
    -Cmh3(ix2,ix3)/(dt*dx3iall(ix3)*dx3all(ix3))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3)    !BC's and pol. time deriv.

    coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)*dx2iall(ix2+1)) ) &
    -Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2)*dx2iall(ix2-1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)*dx3all(ix3+1)) ) &
    -Cm(ix2,ix3-1)*v3(ix2,ix3-1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3)*dx3iall(ix3-1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    ient=ient+1


    !!!ix2+1,ix3 grid point
    ir(ient)=iPhi
    ic(ient)=iPhi+1

    M(ient)=SigPh2(ix2+1,ix3)/(dx2iall(ix2)*dx2all(ix2+1))+gradSigH3(ix2,ix3)/(dx2all(ix2)+dx2all(ix2+1))    !static

    coeff=Cmh2(ix2+1,ix3)/(dt*dx2iall(ix2)*dx2all(ix2+1))
    M(ient)=M(ient)+coeff    !pol. time deriv. terms
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2+1,ix3)    !BC's and pol. time deriv.

    coeff=-Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+2)*dx2iall(ix2+1)) ) &
    -Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)*dx2iall(ix2+1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=-Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) ) &
    -Cm(ix2,ix3-1)*v2(ix2,ix3-1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

    ient=ient+1


    !ix2+2,ix3 grid point
    coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+2)*dx2iall(ix2+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+2

      M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

      ient=ient+1
    end if


    !ix2-2,ix3+1 grid point
    coeff=Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(ix3+1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+lx2-2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !!!ix2,ix3+1 grid point
    ir(ient)=iPhi
    ic(ient)=iPhi+lx2

    M(ient)=SigPh3(ix2,ix3+1)/(dx3iall(ix3)*dx3all(ix3+1))-gradSigH2(ix2,ix3)/(dx3all(ix3)+dx3all(ix3+1))    !static

    coeff=Cmh3(ix2,ix3+1)/(dt*dx3iall(ix3)*dx3all(ix3+1))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3+1)    !BC's and pol. time deriv.

    coeff=-Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+2)*dx3iall(ix3+1)) ) &
    -Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)*dx3iall(ix3+1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    coeff=-Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) ) &
    -Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

    ient=ient+1


    !ix2+2,ix3+1 grid point
    coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3+1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+lx2+2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2-1,ix3+2 grid point
    coeff=-Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    ( (dx3all(ix3+1)+dx3all(ix3+2))*(dx3all(ix3)+dx3all(ix3+1))*(dx2all(ix2)+dx2all(ix2+1)) )
    if (ix3==lx3-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2-1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+2*lx2-1

      M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

      ient=ient+1
    end if


    !ix2,ix3+2 grid point
    coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+2)*dx3iall(ix3+1)) )
    if (ix3==lx3-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+2*lx2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2+1,ix3+2 grid point
    coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) )
    if (ix3==lx3-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx3(ix2+1)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+2*lx2+1

      M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

      ient=ient+1
    end if

  end do loopx2
end do loopx3
!end if


!FIRE UP MUMPS
!if (myid == 0) then
if (debug) print *, 'Filled ',ient-1,' matrix entries.  Initializing MUMPS...'
!end if
mumps_par%COMM = MPI_COMM_WORLD%mpi_val
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)


!LOAD OUR PROBLEM
!if ( myid==0 ) then
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

if (perflag .and. it/=1) then       !used cached permutation
  allocate(mumps_par%PERM_IN(mumps_par%N))
  mumps_par%PERM_IN=mumps_perm
  mumps_par%ICNTL(7)=1
end if

!may solve some memory allocation issues, uncomment if MUMPS throws errors
!about not having enough memory
!mumps_par%ICNTL(14)=50
!end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic2D_polarization')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
!if ( myid==0 ) then
if (debug) print *, 'Now organizing results...'

if (perflag .and. it==1) then
  allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
  mumps_perm=mumps_par%SYM_PERM
end if

elliptic2D_polarization=reshape(mumps_par%RHS,[lx2,lx3])

if (debug) print *, 'Now attempting deallocations...'

deallocate( mumps_par%IRN )
deallocate( mumps_par%JCN )
deallocate( mumps_par%A   )
deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end procedure elliptic2D_polarization


module procedure elliptic2D_polarization_periodic

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS
!-------INCLUDES FULL OF POLARIZATION CURRENT, INCLUDING CONVECTIVE
!-------TERMS.  VELOCITIES SHOULD BE TRIMMED (WITHOUT GHOST CELLS).
!-------THIS VERSION OF THE *INTEGRATED* POTENTIAL SOLVER OBVIATES
!-------ALL OTHERS SINCE A PURELY ELECTRSTATIC FORM CAN BE RECOVERED
!-------BY ZEROING OUT THE INERTIAL CAPACITANCE.
!-------
!-------THIS FORM IS INTENDED TO  WORK WITH CARTESIAN MESHES ONLY.
!-------NOTE THAT THE FULL GRID VARIABLES (X%DX3ALL, ETC.) MUST
!-------BE USED HERE!!!
!-------
!-------THIS FUNCTION WORKS ON A PERIODIC MESH BY USING A CIRCULANT MATRIX
!-------The equation solved by this subroutine is:
!-------
!-------  d/dx2(A dV/dx2) + d/dx3(A dV/dx3) + B dV/dx2 + C dV/x3 + ...
!-------  d/dx2(D d/dt(dV/dx2)) + d/dx3(D d/dt(dV/dx3)) + ...
!-------  d/dx2( D*v2 d^V/dx2^2 + D*v3*d^2V/dx3/dx2 ) + ...
!-------  d/dx3( D*v2 d^2V/dx2/dx3 + D*v3 d^2V/dx3^2 ) = srcterm
!------------------------------------------------------------

real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh2
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: SigPh3
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: Cmh2
real(wp), dimension(1:size(SigP,1),1:size(SigP,2)) :: Cmh3

real(wp) :: coeff    !coefficient for calculating polarization terms
integer :: ix2,ix3,lx2,lx3    !this overwrites the values stored in the grid module, which is fine, but perhaps redundant
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b

type(MUMPS_STRUC) :: mumps_par

integer :: lcount,ix2tmp,ix3tmp

real(wp), dimension(size(SigP,1),size(SigP,2)) :: tmpresults


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then
lx2=size(SigP,1)    !these are full-grid sizes since grid module globals are not in scope
lx3=size(SigP,2)
lPhi=lx2*lx3
!      lent=5*(lx2-2)*(lx3-2)+2*lx2+2*(lx3-2)
!! static model
!      lent=17*(lx2-2)*(lx3A-2)+2*lx2+2*(lx3-2)-3*2*(lx2-2)-3*2*(lx3-2)
!! interior+boundary-x3_adj-x2_adj.  Note that are 3 sets of entries for each adjacent point.
!! The x3 adjacent points do not need to be removed in the case of periodic boundary conditions.
!! This is technicall correct, but a bit misleading, I think.
!! Shouldn't it be lent=17*(lx2-2)*(lx3-2)+2*(lx2-2)+2*lx3-3*2*(lx2-2)-3*2*(lx3-2)???????
!      lent=17*(lx2-2)*(lx3+1-2)+2*(lx3+1)+3*(lx2-2)-3*2*(lx3+1-2)
!! true interior with x3 boundaries which are not treated as interior in periodic solves
!! + add x2 boundaries (note that these are now size lx3+1) + 3 entries for each x3 ghost cell that we are adding
!! - x2_adj (x2 is not periodici, two sets of three points each, note again the larger x3 size as compared to aperiodic solutions).
lent=17*(lx2-2)*(lx3-2)+2*(lx3)+17*2*(lx2-2)-3*2*lx3
!! true interior with x3 boundaries which are not treated as interior in periodic solves
!! + add x2 boundaries (note that these are now size lx3+1) + x3 edge cells (treated here as interior)
!! - x2_adj (x2 is not periodici, two sets of three points each, note again the larger x3 size as compared to aperiodic solutions).
allocate(ir(lent),ic(lent),M(lent),b(lPhi))

if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx2,lx3
if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


!NOTE THAT THESE NEED TO BE PERIODIC IN X3
SigPh2(1,:)= 0
SigPh2(2:lx2,:)=0.5_wp*(SigP(1:lx2-1,:)+SigP(2:lx2,:))
SigPh3(:,1)=0.5_wp*(SigP(:,lx3)+SigP(:,1))    !needs to be left interface value so average of first and last grid point
SigPh3(:,2:lx3)=0.5_wp*(SigP(:,1:lx3-1)+SigP(:,2:lx3))
Cmh2(1,:)= 0
Cmh2(2:lx2,:)=0.5_wp*(Cm(1:lx2-1,:)+Cm(2:lx2,:))
Cmh3(:,1)=0.5_wp*(Cm(:,lx3)+Cm(:,1))
Cmh3(:,2:lx3)=0.5_wp*(Cm(:,1:lx3-1)+Cm(:,2:lx3))


!------------------------------------------------------------
!-------DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED
!-------ASSEMBLED MATRIX INPUT, SEE SECTION 4.5 OF MUMPS USER
!-------GUIDE).
!------------------------------------------------------------
if (debug) print *, 'Loading up matrix entries...'

!LOAD UP MATRIX ELEMENTS
lcount=0
M(:)= 0
b=pack(srcterm,.true.)           !boundaries overwritten later, polarization terms also added later.
ient=1
do ix3=1,lx3
!! note that we have one extra ghost cell now to deal with due to the way we've chosen to implement periodic boundary conditions
  do ix2=1,lx2
    iPhi=lx2*(ix3-1)+ix2
    !! linear index referencing Phi(ix2,ix3) as a column vector.  Also row of big matrix

    if (ix2==1) then
      !! BOTTOM GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vminx2(ix3)
      ient=ient+1
      cycle
    elseif (ix2==lx2) then
      !! TOP GRID POINTS + CORNER
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=1
      b(iPhi)=Vmaxx2(ix3)
      ient=ient+1
      cycle
    endif

    !! TREAT AS AN INTERIOR LOCATION, THIS INCLUDE X3 EDGES NOW SINCE PERIODIC,CIRCULANT
    !! ZZZ - NEED TO WRAP INDICES AROUND:  X 1)
    !! matrix row/column entries; 2) references to dx3i*(anything but ix3);  3) references to conductances/bcs/etc.

    !ix2-1,ix3-2 grid point
    coeff=-Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    ir(ient)=iPhi
    ix3tmp=mod(ix3-2-1+lx3,lx3)+1
    ix2tmp=ix2-1
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    !              ic(ient)=iPhi-2*lx2-1+lPhi-1
    !! add the grid size to wrap the index around (would be negative otherwise),
    !! -1 at end because the last grid opint is actually the same as the first for our implementation
    !            else
    !              ic(ient)=iPhi-2*lx2-1
    !            end if
    M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
    ient=ient+1


    !ix2,ix3-2 grid point
    coeff=-Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)* &
    dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
    ir(ient)=iPhi
    ix3tmp=mod(ix3-2-1+lx3,lx3)+1
    ix2tmp=ix2
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    !              ic(ient)=iPhi-2*lx2+lPhi-1    !again add grid size to wrap to end
    !            else
    !              ic(ient)=iPhi-2*lx2
    !            end if
    M(ient)=coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term
    ient=ient+1


    !ix2+1,ix3-2 grid point
    coeff=Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    ir(ient)=iPhi
    ix3tmp=mod(ix3-2-1+lx3,lx3)+1
    ix2tmp=ix2+1
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
    ient=ient+1


    !ix2-2,ix3-1
    coeff=-Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(mod(ix3-1-1+lx3,lx3)+1)
    else
      ir(ient)=iPhi
      ix3tmp=mod(ix3-1-1+lx3,lx3)+1
      ix2tmp=ix2-2
      ic(ient)=lx2*(ix3tmp-1)+ix2tmp
      !              ic(ient)=iPhi-lx2-2

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2,ix3-1 grid point in ix2,ix3 equation
    ir(ient)=iPhi
    !            ic(ient)=iPhi-lx2
    ix3tmp=mod(ix3-1-1+lx3,lx3)+1
    ix2tmp=ix2
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp

    M(ient)=SigPh3(ix2,ix3)/(dx3iall(ix3)*dx3all(ix3))+gradSigH2(ix2,ix3)/(dx3all(ix3)+dx3all(ix3+1))    !static terms

    coeff=Cmh3(ix2,ix3)/(dt*dx3iall(ix3)*dx3all(ix3))
    M(ient)=M(ient)+coeff   !polarization time derivative terms
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,mod(ix3-1-1+lx3,lx3)+1)
    !! add in polarziation terms that include previous time step potential at this grid point

    coeff=Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3)*dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )+ &
    Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)*dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
    M(ient)=M(ient)+coeff   !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) )+ &
    Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    M(ient)=M(ient)+coeff   !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

    ient=ient+1


    !ix2+2,ix3-1 grid point
    coeff=-Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2+1)+dx2all(ix2+2))*(dx2all(ix2)+dx2all(ix2+1))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(mod(ix3-1-1+lx3,lx3)+1)
    else
      ir(ient)=iPhi
      !              ic(ient)=iPhi-lx2+2
      ix3tmp=mod(ix3-1-1+lx3,lx3)+1
      ix2tmp=ix2+2
      ic(ient)=lx2*(ix3tmp-1)+ix2tmp

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2-2,ix3 grid point
    coeff=-Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)*dx2iall(ix2-1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(ix3)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi-2

      M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

      ient=ient+1
    end if


    !ix2-1,ix3 grid point
    ir(ient)=iPhi
    ic(ient)=iPhi-1

    M(ient)=SigPh2(ix2,ix3)/(dx2iall(ix2)*dx2all(ix2))-gradSigH3(ix2,ix3)/(dx2all(ix2)+dx2all(ix2+1))    !static

    coeff=Cmh2(ix2,ix3)/(dt*dx2iall(ix2)*dx2all(ix2))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2-1,ix3)    !BC's and pol. time deriv.

    coeff=Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2)*dx2iall(ix2-1)) )+ &
    Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)*dx2iall(ix2-1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) )+ &
    Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    M(ient)=M(ient)+coeff     !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

    ient=ient+1


    !ix2,ix3 grid point (main diagonal)
    ir(ient)=iPhi
    ic(ient)=iPhi

    M(ient)=-SigPh2(ix2+1,ix3)/(dx2iall(ix2)*dx2all(ix2+1)) &
    -SigPh2(ix2,ix3)/(dx2iall(ix2)*dx2all(ix2)) &
    -SigPh3(ix2,mod(ix3+1-1,lx3)+1)/(dx3iall(ix3)*dx3all(ix3+1)) &
    -SigPh3(ix2,ix3)/(dx3iall(ix3)*dx3all(ix3))    !static

    coeff=-Cmh2(ix2+1,ix3)/(dt*dx2iall(ix2)*dx2all(ix2+1)) &
    -Cmh2(ix2,ix3)/(dt*dx2iall(ix2)*dx2all(ix2)) &
    -Cmh3(ix2,mod(ix3+1-1,lx3)+1)/(dt*dx3iall(ix3)*dx3all(ix3+1)) &
    -Cmh3(ix2,ix3)/(dt*dx3iall(ix3)*dx3all(ix3))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,ix3)    !BC's and pol. time deriv.

    coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)*dx2iall(ix2+1)) ) &
    -Cm(ix2-1,ix3)*v2(ix2-1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2)*dx2iall(ix2-1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)*dx3all(ix3+1)) ) &
    -Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v3(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3)*dx3iall(mod(ix3-1-1+lx3,lx3)+1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    ient=ient+1


    !ix2+1,ix3 grid point
    ir(ient)=iPhi
    ic(ient)=iPhi+1

    M(ient)=SigPh2(ix2+1,ix3)/(dx2iall(ix2)*dx2all(ix2+1))+gradSigH3(ix2,ix3)/(dx2all(ix2)+dx2all(ix2+1))    !static

    coeff=Cmh2(ix2+1,ix3)/(dt*dx2iall(ix2)*dx2all(ix2+1))
    M(ient)=M(ient)+coeff    !pol. time deriv. terms
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2+1,ix3)    !BC's and pol. time deriv.

    coeff=-Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+2)*dx2iall(ix2+1)) ) &
    -Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)*dx2iall(ix2+1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

    coeff=-Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) ) &
    -Cm(ix2,mod(ix3-1-1+lx3,lx3)+1)*v2(ix2,mod(ix3-1-1+lx3,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3-1)+dx3all(ix3))*(dx2all(ix2)+dx2all(ix2+1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )

    ient=ient+1


    !ix2+2,ix3 grid point
    coeff=Cm(ix2+1,ix3)*v2(ix2+1,ix3)/( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+2)*dx2iall(ix2+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(ix3)
    else
      ir(ient)=iPhi
      ic(ient)=iPhi+2

      M(ient)=coeff    !d/dx2( Cm*v2*d^2/dx2^2(Phi) ) term

      ient=ient+1
    end if


    !ix2-2,ix3+1 grid point
    coeff=Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==2) then
      b(iPhi)=b(iPhi)-coeff*Vminx2(mod(ix3+1-1,lx3)+1)
    else
      ir(ient)=iPhi
      !              ic(ient)=iPhi+lx2-2
      ix3tmp=mod(ix3+1-1,lx3)+1
      ix2tmp=ix2-2
      ic(ient)=lx2*(ix3tmp-1)+ix2tmp

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2,ix3+1 grid point
    ir(ient)=iPhi
    !            ic(ient)=iPhi+lx2
    ix3tmp=mod(ix3+1-1,lx3)+1
    ix2tmp=ix2
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp

    M(ient)=SigPh3(ix2,mod(ix3+1-1,lx3)+1)/ &
    (dx3iall(ix3)*dx3all(ix3+1))-gradSigH2(ix2,ix3)/(dx3all(ix3)+dx3all(ix3+1))    !static

    coeff=Cmh3(ix2,mod(ix3+1-1,lx3)+1)/(dt*dx3iall(ix3)*dx3all(ix3+1))
    M(ient)=M(ient)+coeff    !pol. time deriv.
    b(iPhi)=b(iPhi)+coeff*Phi0(ix2,mod(ix3+1-1,lx3)+1)    !BC's and pol. time deriv.

    coeff=-Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+2)*dx3iall(mod(ix3+1-1,lx3)+1)) ) &
    -Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)*dx3iall(mod(ix3+1-1,lx3)+1)) )
    M(ient)=M(ient)+coeff    !d/dx3( Cm*v3*d^2/dx3^2(Phi) ) term

    coeff=-Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) ) &
    -Cm(ix2-1,ix3)*v3(ix2-1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2-1)+dx2all(ix2))*(dx3all(ix3)+dx3all(ix3+1)) )
    M(ient)=M(ient)+coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

    ient=ient+1


    !ix2+2,ix3+1 grid point
    coeff=Cm(ix2+1,ix3)*v3(ix2+1,ix3)/ &
    ( (dx2all(ix2)+dx2all(ix2+1))*(dx2all(ix2+1)+dx2all(ix2+2))*(dx3all(ix3)+dx3all(ix3+1)) )
    if (ix2==lx2-1) then
      b(iPhi)=b(iPhi)-coeff*Vmaxx2(mod(ix3+1-1,lx3)+1)
    else
      ir(ient)=iPhi
      !              ic(ient)=iPhi+lx2+2
      ix3tmp=mod(ix3+1-1,lx3)+1
      ix2tmp=ix2+2
      ic(ient)=lx2*(ix3tmp-1)+ix2tmp

      M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )

      ient=ient+1
    end if


    !ix2-1,ix3+2 grid point
    !            coeff=-Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    !                  ( (dx3all(ix3+1)+dx3all(ix3+2))*(dx3all(ix3)+dx3all(ix3+1))*(dx2all(ix2)+dx2all(ix2+1)) )
    coeff=-Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &      !mods to wrap indices around
    ( (dx3all(ix3+1)+dx3all(ix3+2))*(dx3all(ix3)+dx3all(ix3+1))*(dx2all(ix2)+dx2all(ix2+1)) )
    ir(ient)=iPhi
    !            if (ix3>=lx3-1) then    !this needs to also handle the case where ix3=lx3!!!  Likewise for statements that follow...
    ix3tmp=mod(ix3+2-1,lx3)+1
    ix2tmp=ix2-1
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    !              ic(ient)=iPhi+2*lx2-1-lPhi+1    !wrap to beginning
    !            else
    !              ic(ient)=iPhi+2*lx2-1
    !            end if
    M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
    ient=ient+1


    !ix2,ix3+2 grid point
    !            coeff=Cm(ix2,ix3+1)*v3(ix2,ix3+1)/ &
    !                  ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+2)*dx3iall(ix3+1)) )
    coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v3(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+2)*dx3iall(mod(ix3+1-1,lx3)+1)) )
    ir(ient)=iPhi
    !            if (ix3>=lx3-1) then
    ix3tmp=mod(ix3+2-1,lx3)+1
    ix2tmp=ix2
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    !              ic(ient)=iPhi+2*lx2-lPhi+1    !subtract grid size to wrap around to the beginning
    !            else
    !              ic(ient)=iPhi+2*lx2
    !            end if
    M(ient)=coeff    !d/dx2( Cm*v3*d^2/dx2dx3(Phi) )
    ient=ient+1


    !ix2+1,ix3+2 grid point
    !            coeff=Cm(ix2,ix3+1)*v2(ix2,ix3+1)/ &
    !                  ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) )
    coeff=Cm(ix2,mod(ix3+1-1,lx3)+1)*v2(ix2,mod(ix3+1-1,lx3)+1)/ &
    ( (dx3all(ix3)+dx3all(ix3+1))*(dx3all(ix3+1)+dx3all(ix3+2))*(dx2all(ix2)+dx2all(ix2+1)) )
    ir(ient)=iPhi
    !            if (ix3>=lx3-1) then     !this should actually work for any case...
    ix3tmp=mod(ix3+2-1,lx3)+1
    ix2tmp=ix2+1
    ic(ient)=lx2*(ix3tmp-1)+ix2tmp
    !              ic(ient)=iPhi+2*lx2+1-lPhi+1
    !            else
    !              ic(ient)=iPhi+2*lx2+1
    !            end if
    M(ient)=coeff    !d/dx3( Cm*v2*d^2/dx3dx2(Phi) )
    ient=ient+1
  end do
end do
!end if


!FIRE UP MUMPS
!if (myid == 0) then
if (debug) print *,  'Debug count:  ',lcount
if (debug) print *, 'Filled ',ient-1,' out of ',lent,' matrix entries for solving ',iPhi,' of ',lPhi, &
' unknowns.  Initializing MUMPS...'
if (ient-1 /= lent) error stop 'Incorrect number of matrix entries filled in potential solve!!!'

!end if
mumps_par%COMM = MPI_COMM_WORLD%mpi_val
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)


!LOAD OUR PROBLEM
!if ( myid==0 ) then
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

if (perflag .and. it/=1) then       !used cached permutation
  allocate(mumps_par%PERM_IN(mumps_par%N))
  mumps_par%PERM_IN=mumps_perm
  mumps_par%ICNTL(7)=1
end if

!mumps_par%ICNTL(14)=50
!end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic2D_polarization_periodic')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
!if ( myid==0 ) then
if (debug) print *, 'Now organizing results...'

if (perflag .and. it==1) then
  allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
  mumps_perm=mumps_par%SYM_PERM
end if

!IF WE HAVE DONE A PERIODIC SOLVE, THE LAST GRID POINT NEEDS TO BE IGNORED WHEN WE RESHAPE THE POTENTIAL ARRAY.

tmpresults=reshape(mumps_par%RHS,[lx2,lx3])
elliptic2D_polarization_periodic=tmpresults(1:lx2,1:lx3)    !sort of superfluous now that the solve size is the same as the grid

if (debug) print *, 'Now attempting deallocations...'

deallocate( mumps_par%IRN )
deallocate( mumps_par%JCN )
deallocate( mumps_par%A   )
deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end procedure elliptic2D_polarization_periodic


! FIXME:  x3 below really refers to whatever the second non-singleton dimension is...
module procedure elliptic2D_cart
!! SOLVE IONOSPHERIC POTENTIAL EQUATION IN 2D USING MUMPS.
!!
!! ASSUME THAT:
!! * WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD LINE
!! * POTENTIAL VARIES IN X1 AND X3. X2 IS NOMINALLY JUST ONE ELEMENT.
!! * LEFT AND RIGHT BOUNDARIES (IN X3) USE DIRICHLET BOUNDARY CONDITIONS
!! * TOP (ALTITUDE) CAN BE NEUMANN OR DIRICHLET.
!! * BOTTOM (ALTITUDE) IS ALWAYS Neumann zero current
!!
!! This subroutine solves equations of the form:
!!
!!    d/dx1(sig0 dV/dx1) + d/dx3(sigP dV/dx3) = srcterm
!!
!! The boundary conditions arrays provided to this procedure are assumed to be
!! in units of Volts (dirichlet) or V/m (Neumann), meaning any currents need
!! to be converted into potential normal derivatives prior to calling this.


real(wp), dimension(:,:), allocatable :: sig0h1
real(wp), dimension(:,:), allocatable :: sigPh3

integer :: ix1,ix3,lx1,lx3
integer :: lx2,l2nddim
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b

integer :: ibnd,ldirichx1,lneux1,ldirichx3,lneux3
logical :: flag2
type(MUMPS_STRUC) :: mumps_par


!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then


! system size
lx1=size(sig0,1)
lx2=size(sig0,2)
lx3=size(sig0,3)

if (lx2==1) then
  flag2=.false.
  l2nddim=lx3
else if (lx3==1) then
  flag2=.true.
  l2nddim=lx2
else
  error stop ' elliptic2d_cart -> could not determine which dimensions of arrays to use!!!'
end if

!number of neumann boundaries (x1 and x3) used in this solve
lneux1=0
do ibnd=1,2
  if (flagsdirich(ibnd)==0) lneux1=lneux1+1
end do
ldirichx1=2-lneux1
lneux3=0
do ibnd=3,4
  if (flagsdirich(ibnd)==0) lneux3=lneux3+1
end do
ldirichx3=2-lneux3

! count the number of matrix entries we need to fill (will depend on type of boundary conditions chosen
lPhi=lx1*l2nddim
!if (flagsdirich==0) then
!  lent=5*(lx1-2)*(lx3-2)+2*lx1+2*(lx3-2)+2*lx3    !first +1 for Neumann bottom, second for Neumann top
!else
!  !        lent=5*(lx1-2)*(lx3-2)+2*lx1+2*(lx3-2)+1    !first +1 for Neumann bottom
!  lent=5*(lx1-2)*(lx3-2)+2*(lx1-2)+2*lx3+lx3
!end if

! count matrix entries as follows:  interior points (5 entries each) + # dirich x1 * size + # neumann x1 * size + # dirich x3 * size sans corners + # neumann * size sans corners
lent=5*(lx1-2)*(l2nddim-2) + ldirichx1*l2nddim + lneux1*2*l2nddim + ldirichx3*(lx1-2) + lneux3*2*(lx1-2)


! allocate space for our problem
allocate(ir(lent),ic(lent),M(lent),b(lPhi))
if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx1,l2nddim
if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent
if (debug) print*, 'Number of Neumann boundaries:  ', lneux1,lneux3
if (debug) print*, 'Number of Dirichlet boundaries:  ', ldirichx1,ldirichx3


! conductivities need to be averaged to the cell interfaces for the FDE we use
allocate(sig0h1(lx1,l2nddim),sigPh3(lx1,l2nddim))
sig0h1(1,:)=0
sigPh3(:,1)=0
if (flag2) then
  sigPh3(:,2:l2nddim)=0.5_wp*(sigP(:,1:l2nddim-1,1)+sigP(:,2:l2nddim,1))
  sig0h1(2:lx1,:)=0.5_wp*(sig0(1:lx1-1,:,1)+sig0(2:lx1,:,1))
else
  sigPh3(:,2:l2nddim)=0.5_wp*(sigP(:,1,1:l2nddim-1)+sigP(:,1,2:l2nddim))
  sig0h1(2:lx1,:)=0.5_wp*(sig0(1:lx1-1,1,:)+sig0(2:lx1,1,:))
end if

! fill elements of matrix to be solved.  we use centralized assembled matrix input as described in mumps user manual section 4.5
! all of the logic of inverted vs. noninverted grids has been exported to the parent routine; leaving this as a pure applied
! math procedure with no specific knowledgeo of the ionospheric problem.
M(:)=0
b=pack(srcterm,.true.)           !boundaries overwritten later
ient=1
do ix3=1,l2nddim
  do ix1=1,lx1
    iPhi=lx1*(ix3-1)+ix1     !linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

    if (ix1==1) then    !! (LOGICAL) BOTTOM GRID POINTS
      if (flagsdirich(1)/=0) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        if (flag2) then
          b(iPhi)=Vminx1(ix3,1)  ! new routines always map min/max as specified in input files
        else
          b(iPhi)=Vminx1(1,ix3)  ! new routines always map min/max as specified in input files
        end if
        ient=ient+1
      else
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1/dx1(2)
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=1/dx1(2)
        if (flag2) then
          b(iPhi)=Vminx1(ix3,1)
        else
          b(iPhi)=Vminx1(1,ix3)
        end if
        ient=ient+1
      end if
    elseif (ix1==lx1) then    !(LOGICAL) TOP GRID POINTS + CORNER
      if (flagsdirich(2)/=0) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        if (flag2) then
          b(iPhi)=Vmaxx1(ix3,1)
        else
          b(iPhi)=Vmaxx1(1,ix3)
        end if
        ient=ient+1
      else
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=-1/dx1(lx1)
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1/dx1(lx1)
        if (flag2) then
          b(iPhi)=Vmaxx1(ix3,1)
        else
          b(iPhi)=Vmaxx1(1,ix3)
        end if
        ient=ient+1
      end if
    elseif (ix3==1) then      !LEFT BOUNDARY
      if (flagsdirich(3)/=0) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vminx3(ix1,1)
        ient=ient+1
      else
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1/dx3all(2)
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=1/dx3all(2)
        b(iPhi)=Vminx3(ix1,1)
        ient=ient+1
      end if
    elseif (ix3==l2nddim) then    !RIGHT BOUNDARY
      if (flagsdirich(4)/=0) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1.0
        b(iPhi)=Vmaxx3(ix1,1)
        ient=ient+1
      else
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=-1/dx3all(l2nddim)
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1/dx3all(l2nddim)
        b(iPhi)=Vmaxx3(ix1,1)
        ient=ient+1
      end if
    else                      !INTERIOR
      !ix1,ix3-1 grid point in ix1,ix3 equation
      ir(ient)=iPhi
      ic(ient)=iPhi-lx1
      M(ient)=sigPh3(ix1,ix3)/(dx3iall(ix3)*dx3all(ix3))
      ient=ient+1

      !ix1-1,ix3 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi-1
      M(ient)=sig0h1(ix1,ix3)/(dx1i(ix1)*dx1(ix1))
      ient=ient+1

      !ix1,ix3 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi
      M(ient)=-sig0h1(ix1+1,ix3)/(dx1i(ix1)*dx1(ix1+1)) &
      -sig0h1(ix1,ix3)/(dx1i(ix1)*dx1(ix1)) &
      -sigPh3(ix1,ix3+1)/(dx3iall(ix3)*dx3all(ix3+1)) &
      -sigPh3(ix1,ix3)/(dx3iall(ix3)*dx3all(ix3))
      ient=ient+1

      !ix1+1,ix3 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi+1
      M(ient)=sig0h1(ix1+1,ix3)/(dx1i(ix1)*dx1(ix1+1))
      ient=ient+1

      !ix1,ix3+1 grid point
      ir(ient)=iPhi
      ic(ient)=iPhi+lx1
      M(ient)=sigPh3(ix1,ix3+1)/(dx3iall(ix3)*dx3all(ix3+1))
      ient=ient+1
    end if
  end do
end do
!end if
if (debug) print *, 'Number of entries used:  ',ient-1


!FIRE UP MUMPS
mumps_par%COMM = MPI_COMM_WORLD%mpi_val
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)


!LOAD OUR PROBLEM
!if ( myid==0 ) then
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

if (perflag .and. it/=1) then       !used cached permutation
  if (debug) print *, 'Using a previously stored permutation'
  allocate(mumps_par%PERM_IN(mumps_par%N))
  mumps_par%PERM_IN = mumps_perm
  mumps_par%ICNTL(7) = 1
end if

!may solve some memory allocation issues, uncomment if MUMPS throws errors
!about not having enough memory
! e.g. INFO(1) = -9
! however this error may also mean there is a deeper problem with the code.
!,mumps_par%ICNTL(14) = 50
!end if


!> SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic2D_cart')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors).  WOW - this halves execution
! time on some big 2048*2048 solves!!!
!if ( myid==0 ) then
if (debug) print *, 'Now organizing results...'

if (perflag .and. it==1) then
  if (debug) print *, 'Storing ordering for future time step use...'
  allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
  mumps_perm=mumps_par%SYM_PERM
end if

if (flag2) then
  elliptic2D_cart=reshape(mumps_par%RHS,[lx1,l2nddim,1])
else
  elliptic2D_cart=reshape(mumps_par%RHS,[lx1,1,l2nddim])
end if

if (debug) print *, 'Now attempting deallocations...'

deallocate( mumps_par%IRN )
deallocate( mumps_par%JCN )
deallocate( mumps_par%A   )
deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

deallocate(sig0h1,sigPh3)

end procedure elliptic2D_cart


end submodule elliptic2d
