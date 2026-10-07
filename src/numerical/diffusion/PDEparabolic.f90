module PDEparabolic

!! a module for use for solving parabolic partial differential equations


!> banded and tridiagonal solvers, for now we just take everything to be banded...
use phys_consts, only: wp
use vendor_lapack95, only: gbsv!,gtsv
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit

implicit none (type, external)

private
public :: TRBDF21D, backEuler1D


contains


function TRBDF21D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,BCtype,dx1,dx1i,Tsmin_mid,Tsmax_mid,increments)

!! SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT
!! GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!! THEY ARE PASSED INTO THIS ROUTINE.
!!
!! FORM OF EQUATION SOLVED IS:
!!  dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
!!
!!  NOTE: UPON FURTHER REVIEW I THINK THE FORM SOLVED IS ACTUALLY:
!!  dT/dt = A T + B dT/dx + C d/dx(D dT/dx) + E
!!
!! We assume the user is providing the type of boundary conditions and that
!! any Neumann boundary conditions are interpreted as diffs in state variable
!! to be solved (this avoids needing solvers to have to do additional calculations
!! using auxiliary variables not in scope.

real(wp), dimension(:), intent(in) :: A,B,C,D,E
real(wp), dimension(:), intent(in) :: Ts
real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
! Existing callers prescribe stage-constant boundaries. Time-varying Dirichlet
! or Neumann data must also supply the value at t+dt/2 for the TR stage.
real(wp), intent(in), optional :: Tsmin_mid,Tsmax_mid
! Integrated reaction, drift, conduction, explicit source and imposed endpoint
! increments. These are evaluated from the numerical stage states, not from
! the interior residual. Endpoints are algebraic boundary reservoirs.
real(wp), intent(out), optional :: increments(size(Ts),5)
integer, dimension(2), intent(in) :: BCtype  !=0 dirichlet; =1 neumann
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
integer, parameter :: ll=2                   !number of lower diagonals

real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh
integer :: ix1,lx1,info

real(wp), dimension(size(Ts)) :: TR

real(wp), dimension(size(Ts)) :: TRBDF21D

!> ORGANIZE SIZES AND THERMAL CONDUCTIVITY
lx1=size(Ts)
Dh(1)=0
Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
!    TR(:)=Ts(:)/dt+E(:)
!! boundaries to be overwritten later...  This is now done for each grid point in a separate statement


!! ## TR HALF STEP:  DEFINE A MATRIX USING BANDED STORAGE


! ZZZ - check whether diriclet or neumann...
!> MINX1 BOUNDARY (DIRICHLET)
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=0
  M(ll+1,ix1+2)=0
  TR(ix1)=Tsminx1
else
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=-1
  M(ll+1,ix1+2)=0
  TR(ix1)=Tsminx1 !used to be 0, but we change it since we are doing the difference.
end if


!> FIRST INTERIOR GRID POINT
ix1=2

!> ix1-1
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2

!> ix1
M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2

!> ix1+1, super-diag.
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
M(ll+1,ix1+2)=0
TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1) &
  -M(ll+1,ix1+2)*Ts(ix1+2)


!> INTERIOR GRID POINTS
do ix1=3,lx1-2
!! removed do concurrent to avoid oneAPI compiler ICE #5623
!! do concurrent OK because only indexing already defined things
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
  M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &    !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
  TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
    -M(ll+5,ix1-2)*Ts(ix1-2) &
    -M(ll+4,ix1-1)*Ts(ix1-1) &
    -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
    -M(ll+2,ix1+1)*Ts(ix1+1) &
    -M(ll+1,ix1+2)*Ts(ix1+2)
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
M(ll+3,ix1)=1/(dt/2)-A(ix1)/2 &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2 &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))/2
TR(ix1)=Ts(ix1)/(dt/2)+E(ix1) &
  -M(ll+5,ix1-2)*Ts(ix1-2) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1)


! ZZZ - check whether dirichlet or neumann...
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  TR(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
  TR(ix1)=Tsmaxx1
end if


if (present(Tsmin_mid)) TR(1)=Tsmin_mid
if (present(Tsmax_mid)) TR(lx1)=Tsmax_mid

!! ### TR HALF STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TR,kl=2,info=info)
call check_solver_info(info,'TRBDF21D TR stage')



!! ## BDF2 STEP:  DEFINE A MATRIX USING BANDED STORAGE

!ZZZ - check whether D or N
!> MINX1 BOUNDARY (DIRICHLET)
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=0
  M(ll+1,ix1+2)=0
  TRBDF21D(ix1)=Tsminx1
else
  M(ll+3,ix1)=1
  M(ll+2,ix1+1)=-1
  M(ll+1,ix1+2)=0
  TRBDF21D(ix1)=Tsminx1
end if


!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/(dt/3)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2)=0
TRBDF21D(ix1)=E(ix1) &
  -1/3._wp*Ts(ix1)/(dt/3) &
  +4/3._wp*TR(ix1)/(dt/3)


!> INTERIOR GRID POINTS
do ix1=3,lx1-2
!! removed do concurrent to avoid oneAPI compiler ICE #5623
  M(ll+5,ix1-2)=0                                               !ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1)=1/(dt/3)-A(ix1) &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &    !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2)=0                                               !ix1+2 grid point
  TRBDF21D(ix1)=E(ix1) &
    -1/3._wp*Ts(ix1)/(dt/3) &
    +4/3._wp*TR(ix1)/(dt/3)
end do


!LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/(dt/3)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
TRBDF21D(ix1)=E(ix1) &
  -1/3._wp*Ts(ix1)/(dt/3) &
  +4/3._wp*TR(ix1)/(dt/3)

!check whether D or N
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  TRBDF21D(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
  TRBDF21D(ix1)=Tsmaxx1
end if


!! ## BDF2 STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TRBDF21D,kl=2,info=info)
call check_solver_info(info,'TRBDF21D BDF2 stage')

if (present(increments)) then
  ! TR uses dt/2 followed by the BDF formula with dt/3. Eliminating the
  ! intermediate state gives dt/3 times each of L(T0), L(TR), L(Tfinal).
  increments=dt/3*(parabolic_terms(Ts,A,B,C,D,E,dx1,dx1i) &
                 +parabolic_terms(TR,A,B,C,D,E,dx1,dx1i) &
                 +parabolic_terms(TRBDF21D,A,B,C,D,E,dx1,dx1i))
  increments(1,5)=TRBDF21D(1)-Ts(1)
  increments(lx1,5)=TRBDF21D(lx1)-Ts(lx1)
endif

end function TRBDF21D


function backEuler1D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,BCtype,dx1,dx1i,coeffs,rhs,increments)

!------------------------------------------------------------
!-------SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  FORM OF EQUATION
!-------SOLVED IS:
!------- dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
!------------------------------------------------------------

real(wp), dimension(:), intent(in) :: A,B,C,D,E
real(wp), dimension(:), intent(in) :: Ts
real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
integer, dimension(2), intent(in) :: BCtype
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
real(wp), dimension(:,:), intent(inout), optional :: coeffs
!! intent(out)
real(wp), dimension(:), intent(inout), optional :: rhs
!! intent(out)
real(wp), intent(out), optional :: increments(size(Ts),5)

integer, parameter :: ll=2                   !number of lower diagonals
real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh
real(wp), dimension(size(Ts)) :: backEuler1D
integer :: ix1,lx1,info

!------------------------------------------------------------
!-------DEFINE A MATRIX USING BANDED STORAGE
!------------------------------------------------------------
lx1=size(Ts)
Dh(1)=0
Dh(2:lx1)=0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
backEuler1D(:)=Ts(:)/dt+E(:)                !boundaries to be overwritten later...


! check whether D or N
!> MINX1 BOUNDARY, Dirichlet BCS
ix1=1
if (BCtype(1)==0) then
  M(ll+3,ix1)=1       !main diagonal denoted temperature at this grid point... 1*Ts,i=Tsminx1
  M(ll+2,ix1+1)=0     !1st super diagonal
  M(ll+1,ix1+2)=0     !2nd super diagonal
  backEuler1D(ix1)=Tsminx1
else
  M(ll+3,ix1)=1       !main diagonal denoted temperature at this grid point... 1*Ts,i=Tsminx1
  M(ll+2,ix1+1)=-1    !1st super diagonal
  M(ll+1,ix1+2)=0     !2nd super diagonal
!  backEuler1D(ix1)=0
  backEuler1D(ix1)=Tsminx1
end if

!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1, sub-diaginal
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2)=0


!> INTERIOR GRID POINTS
do ix1=3,lx1-2
!! removed do concurrent to avoid oneAPI compiler ICE #5623
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
           -B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2)=0
M(ll+4,ix1-1)=-C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1/dt-A(ix1) &                                     !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -B(ix1)/(dx1(ix1+1)+dx1(ix1))

! check whether D or N
!> MAXX1 BOUNDARY
ix1=lx1
if (BCtype(2)==0) then
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=0
  M(ll+3,ix1)=1
  backEuler1D(ix1)=Tsmaxx1
else
  M(ll+5,ix1-2)=0
  M(ll+4,ix1-1)=-1
  M(ll+3,ix1)=1
!  backEuler1D(ix1)=0
  backEuler1D(ix1)=Tsmaxx1
end if

!> in case we want to output the right-hand side of the system; has to be done here before
!   the output argument (which stores rhs) is overwritten by the solution.
if (present(rhs)) then
  rhs(1:lx1)=backEuler1D(1:lx1)
end if

!! ## DO SOME STUFF TO CALL LAPACK'S BANDED SOLVER
!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,backEuler1D,kl=2,info=info)
call check_solver_info(info,'backEuler1D')

if (present(increments)) then
  increments=dt*parabolic_terms(backEuler1D,A,B,C,D,E,dx1,dx1i)
  increments(1,5)=backEuler1D(1)-Ts(1)
  increments(lx1,5)=backEuler1D(lx1)-Ts(lx1)
endif

!> this is for if one wants to output the matrix bands for testing purposes
if (present(coeffs)) then
  coeffs(1:lx1,1)=M(ll+2,1:lx1)
  coeffs(1:lx1,2)=M(ll+3,1:lx1)
  coeffs(1:lx1,3)=M(ll+4,1:lx1)
end if

end function backEuler1D

subroutine check_solver_info(info,stage)
  integer, intent(in) :: info
  character(*), intent(in) :: stage

  if (info == 0) return
  write(stderr,'(a,a,a,i0)') 'PDEparabolic: ',stage,' gbsv INFO=',info
  if (info < 0) error stop 'PDEparabolic: invalid LAPACK argument'
  error stop 'PDEparabolic: singular diffusion matrix'
end subroutine check_solver_info

pure function parabolic_terms(T,A,B,C,D,E,dx,dxi) result(terms)
  real(wp), intent(in) :: T(:),A(:),B(:),C(:),D(:),E(:),dx(0:),dxi(:)
  real(wp) :: terms(size(T),5),left,right
  integer :: i
  terms=0
  do i=2,size(T)-1
    left=0.5_wp*(D(i-1)+D(i))*(T(i)-T(i-1))/dx(i)
    right=0.5_wp*(D(i)+D(i+1))*(T(i+1)-T(i))/dx(i+1)
    terms(i,1)=A(i)*T(i)
    terms(i,2)=B(i)*(T(i+1)-T(i-1))/(dx(i+1)+dx(i))
    terms(i,3)=C(i)*(right-left)/dxi(i)
    terms(i,4)=E(i)
  enddo
end function parabolic_terms

end module PDEparabolic
