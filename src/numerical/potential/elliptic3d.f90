submodule (PDEelliptic) elliptic3d

implicit none (type, external)

contains

module procedure elliptic3D_cart

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.  This version requires that one compute the
!-------Coefficients of the problem beforehand and solves the equation:
!-------
!-------A d^2V/dx1^2 + B d^2V/dx2^2 + C d^2V/dx3^2 + D dV/dx1 + E dV/dx2 + F dV/dx3 = srcterm
!------------------------------------------------------------

integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b

type (MUMPS_STRUC) :: mumps_par

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then
lx1=size(Ac,1)
lx2=size(Ac,2)
lx3=size(Ac,3)
lPhi=lx1*lx2*lx3

lent=7*(lx1-2)*(lx2-2)*(lx3-2)                                                 !interior entries
lent=lent+2*(lx1-2)*(lx2-2)+2*(lx2-2)*(lx3-2)+2*(lx1-2)*(lx3-2)                !6 faces of cube
lent=lent+4*(lx1-2)+4*(lx2-2)+4*(lx3-2)                                        !12 edges
lent=lent+8                                                                    !8 corners, now we have total nonzero entries
lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom

if (flagdirich==0) lent=lent+lx2*lx3
!! more entries if Neumann on top

allocate(ir(lent),ic(lent),M(lent),b(lPhi))

if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


!! DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED ASSEMBLED MATRIX INPUT
!! SEE SECTION 4.5 OF MUMPS USER GUIDE).

!> LOAD UP MATRIX ELEMENTS
M(:)=0
b=pack(srcterm,.true.)           !boundaries overwritten later
ient=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1
      !! linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

      if (ix1==1) then
        !! BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient) = 1
        !            b(iPhi)=Vminx1(ix3)
        b(iPhi) = 0
        !! force bottom current to zero
        ient = ient+1
      elseif (ix1==lx1) then
        !! TOP GRID POINTS + CORNER
        if (flagdirich/=0) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-1
          M(ient)=-1/dx1(lx1)
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1/dx1(lx1)
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        end if
      elseif (ix2==1) then
        !! LEFT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vminx2(ix1,ix3)
        ient=ient+1
      elseif (ix2==lx2) then
        !! RIGHT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vmaxx2(ix1,ix3)
        ient=ient+1
      elseif (ix3==1) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vminx3(ix1,ix2)
        ient=ient+1
      elseif (ix3==lx3) then
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vmaxx3(ix1,ix2)
        ient=ient+1
      else
        !! INTERIOR
        !> ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1*lx2
        M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3)/dx3iall(ix3)-Ec(ix1,ix2,ix3)/(dx3all(ix3+1)+dx3all(ix3))
        ient=ient+1

        !> ix1,ix2-1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2)/dx2iall(ix2)-Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
        ient=ient+1

        !> ix1-1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1)/dx1i(ix1)-Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
        ient=ient+1

        !> ix1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-Ac(ix1,ix2,ix3)*(1/dx2all(ix2+1)/dx2iall(ix2)+1/dx2all(ix2)/dx2iall(ix2))- &
              Bc(ix1,ix2,ix3)*(1/dx3all(ix3+1)/dx3iall(ix3)+1/dx3all(ix3)/dx3iall(ix3))- &
              Cc(ix1,ix2,ix3)*(1/dx1(ix1+1)/dx1i(ix1)+1/dx1(ix1)/dx1i(ix1))
        ient=ient+1

        !> ix1+1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1+1)/dx1i(ix1)+Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
        ient=ient+1

        !> ix1,ix2+1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2+1)/dx2iall(ix2)+Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
        ient=ient+1

        !> ix1,ix2,ix3+1
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1*lx2
        M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3+1)/dx3iall(ix3)+Ec(ix1,ix2,ix3)/(dx3all(ix3+1)+dx3all(ix3))
        ient=ient+1
      end if
    end do
  end do
end do

if (debug) print *, 'Number of entries used:  ',ient-1


!> INIT MUMPS

mumps_par%COMM = MPI_COMM_WORLD%mpi_val
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)

!LOAD OUR PROBLEM
if (debug) print*, 'Loading mumps problem...'
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

if (debug) print*, 'Dealing with permutation',perflag,it
if (perflag .and. it/=1) then       !used cached permutation, but only gets tested on first time step
  allocate(mumps_par%PERM_IN(mumps_par%N))
  mumps_par%PERM_IN=mumps_perm
  mumps_par%ICNTL(7)=1
end if


if (debug) print*, 'Setting memory relaxation...'
!3D solves very often need better memory relaxation
mumps_par%ICNTL(14)=500
!end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
if (debug) print*, 'Executing solve...'
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic3D_cart')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
!if ( myid==0 ) then
if (debug) print *, 'Now organizing results...'

if (perflag .and. it==1) then
  allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
  mumps_perm=mumps_par%SYM_PERM
end if

elliptic3D_cart=reshape(mumps_par%RHS,[lx1,lx2,lx3])

if (debug) print *, 'Now attempting deallocations...'

deallocate( mumps_par%IRN )
deallocate( mumps_par%JCN )
deallocate( mumps_par%A   )
deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end procedure elliptic3D_cart


module procedure elliptic3D_cart_periodic

!------------------------------------------------------------
!-------SOLVE IONOSPHERIC POTENTIAL EQUATION IN 3D USING MUMPS
!-------ASSUME THAT WE ARE RESOLVING THE POTENTIAL ALONG THE FIELD
!-------LINE.  THIS IS MOSTLY INEFFICIENT/UNWORKABLE FOR MORE THAN 1M
!-------GRID POINTS.  This version requires that one compute the
!-------Coefficients of the problem beforehand and solves the equation:
!-------
!-------A d^2V/dx1^2 + B d^2V/dx2^2 + C d^2V/dx3^2 + D dV/dx1 + E dV/dx2 + F dV/dx3 = srcterm
!------------------------------------------------------------

integer :: ix1,ix2,ix3,lx1,lx2,lx3
integer :: lPhi,lent
integer :: iPhi,ient
integer, dimension(:), allocatable :: ir,ic
real(wp), dimension(:), allocatable :: M
real(wp), dimension(:), allocatable :: b
integer :: ix3prev,ix3next

type (MUMPS_STRUC) :: mumps_par

!ONLY ROOT NEEDS TO ASSEMBLE THE MATRIX
!if (myid==0) then
lx1=size(Ac,1)
lx2=size(Ac,2)
lx3=size(Ac,3)
lPhi=lx1*lx2*lx3

lent=7*(lx1-2)*(lx2-2)*(lx3)                                                   !interior entries include periodic x3 (so all x3 interior)
lent=lent+2*(lx2-2)*(lx3)+2*(lx1-2)*(lx3)                                      !4 faces of periodic cube
lent=lent+4*(lx3)                                                              !4 edges
lent=lent+lx2*lx3                                                              !entries to deal with Neumann conditions on bottom

if (flagdirich==0) lent=lent+lx2*lx3
!! more entries if Neumann on top

allocate(ir(lent),ic(lent),M(lent),b(lPhi))

if (debug) print *, 'MUMPS will attempt a solve of size:  ',lx1,lx2,lx3
if (debug) print *, 'Total unknowns and nonzero entries in matrix:  ',lPhi,lent


!! DEFINE A MATRIX USING SPARSE STORAGE (CENTRALIZED ASSEMBLED MATRIX INPUT
!! SEE SECTION 4.5 OF MUMPS USER GUIDE).

!> LOAD UP MATRIX ELEMENTS
M(:)=0
b=pack(srcterm,.true.)           !boundaries overwritten later
ient=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      iPhi=lx1*lx2*(ix3-1)+lx1*(ix2-1)+ix1
      !! linear index referencing Phi(ix1,ix3) as a column vector.  Also row of big matrix

      if (ix1==1) then
        !! BOTTOM GRID POINTS + CORNER, USE NEUMANN HERE, PRESUMABLY ZERO
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-1
        ient=ient+1
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient) = 1
        !            b(iPhi)=Vminx1(ix3)
        b(iPhi) = 0
        !! force bottom current to zero
        ient = ient+1
      elseif (ix1==lx1) then
        !! TOP GRID POINTS + CORNER
        if (flagdirich/=0) then
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        else
          ir(ient)=iPhi
          ic(ient)=iPhi-1
          M(ient)=-1/dx1(lx1)
          ient=ient+1
          ir(ient)=iPhi
          ic(ient)=iPhi
          M(ient)=1/dx1(lx1)
          b(iPhi)=Vmaxx1(ix2,ix3)
          ient=ient+1
        end if
      elseif (ix2==1) then
        !! LEFT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vminx2(ix1,ix3)
        ient=ient+1
      elseif (ix2==lx2) then
        !! RIGHT BOUNDARY
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=1
        b(iPhi)=Vmaxx2(ix1,ix3)
        ient=ient+1
      else
        ! check if we are going to circulate off the grid
        if (ix3==lx3) then
          ix3next=1
        else
          ix3next=ix3+1
        end if
        if (ix3==1) then
          ix3prev=lx3
        else
          ix3prev=ix3-1
        end if

        !! INTERIOR
        !> ix1,ix2,ix3-1 grid point in ix1,ix2,ix3 equation
        ir(ient)=iPhi
        ic(ient)=lx1*lx2*(ix3prev-1)+lx1*(ix2-1)+ix1
        M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3)/dx3iall(ix3)-Ec(ix1,ix2,ix3)/(dx3all(ix3next)+dx3all(ix3))
        ient=ient+1

        !> ix1,ix2-1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-lx1
        M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2)/dx2iall(ix2)-Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
        ient=ient+1

        !> ix1-1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi-1
        M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1)/dx1i(ix1)-Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
        ient=ient+1

        !> ix1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi
        M(ient)=-Ac(ix1,ix2,ix3)*(1/dx2all(ix2+1)/dx2iall(ix2)+1/dx2all(ix2)/dx2iall(ix2))- &
              Bc(ix1,ix2,ix3)*(1/dx3all(ix3next)/dx3iall(ix3)+1/dx3all(ix3)/dx3iall(ix3))- &
              Cc(ix1,ix2,ix3)*(1/dx1(ix1+1)/dx1i(ix1)+1/dx1(ix1)/dx1i(ix1))
        ient=ient+1

        !> ix1+1,ix2,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+1
        M(ient)=Cc(ix1,ix2,ix3)/dx1(ix1+1)/dx1i(ix1)+Fc(ix1,ix2,ix3)/(dx1(ix1+1)+dx1(ix1))
        ient=ient+1

        !> ix1,ix2+1,ix3
        ir(ient)=iPhi
        ic(ient)=iPhi+lx1
        M(ient)=Ac(ix1,ix2,ix3)/dx2all(ix2+1)/dx2iall(ix2)+Dc(ix1,ix2,ix3)/(dx2all(ix2+1)+dx2all(ix2))
        ient=ient+1

        !> ix1,ix2,ix3+1
        ir(ient)=iPhi
        ic(ient)=lx1*lx2*(ix3next-1)+lx1*(ix2-1)+ix1
        M(ient)=Bc(ix1,ix2,ix3)/dx3all(ix3next)/dx3iall(ix3)+Ec(ix1,ix2,ix3)/(dx3all(ix3next)+dx3all(ix3))
        ient=ient+1
      end if
    end do
  end do
end do

if (debug) print *, 'Number of entries used:  ',ient-1


!> INIT MUMPS
mumps_par%COMM = MPI_COMM_WORLD%mpi_val
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1

call MUMPS_exec(mumps_par)

call quiet_mumps(mumps_par)

!LOAD OUR PROBLEM
if (debug) print*, 'Loading mumps problem...'
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

if (debug) print*, 'Dealing with permutation',perflag,it
if (perflag .and. it/=1) then       !used cached permutation, but only gets tested on first time step
  allocate(mumps_par%PERM_IN(mumps_par%N))
  mumps_par%PERM_IN=mumps_perm
  mumps_par%ICNTL(7)=1
end if


if (debug) print*, 'Setting memory relaxation...'
!3D solves very often need better memory relaxation
mumps_par%ICNTL(14)=500
!end if


!SOLVE (ALL WORKERS NEED TO SEE THIS CALL)
if (debug) print*, 'Executing solve...'
mumps_par%JOB = 6

call MUMPS_exec(mumps_par)

call check_mumps_status(mumps_par, 'elliptic3D_cart')

!STORE PERMUTATION USED, SAVE RESULTS, CLEAN UP MUMPS ARRAYS
!(can save ~25% execution time and improves scaling with openmpi
! ~25% more going from 1-2 processors)
!if ( myid==0 ) then
if (debug) print *, 'Now organizing results...'

if (perflag .and. it==1) then
  allocate(mumps_perm(mumps_par%N))     !we don't have a corresponding deallocate statement
  mumps_perm=mumps_par%SYM_PERM
end if

elliptic3D_cart_periodic=reshape(mumps_par%RHS,[lx1,lx2,lx3])

if (debug) print *, 'Now attempting deallocations...'

deallocate( mumps_par%IRN )
deallocate( mumps_par%JCN )
deallocate( mumps_par%A   )
deallocate( mumps_par%RHS )
!end if

mumps_par%JOB = -2

call MUMPS_exec(mumps_par)

end procedure elliptic3D_cart_periodic

end submodule elliptic3d
