module PDEparabolic
!! a module for use for solving parabolic partial differential equations
!! banded and tridiagonal solvers, for now we just take everything to be banded...

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only: wp
use vendor_lapack95, only: gbsv

implicit none (type, external)

private
public :: TRBDF21D, backEuler1D

contains


function TRBDF21D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)

!! SOLVE A 1D DIFFUSION PROBLEM.  IT IS EXPECTED THAT
!! GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!! THEY ARE PASSED INTO THIS ROUTINE.
!!
!! FORM OF EQUATION SOLVED IS:
!!  dT/dt + A T + B dT/dx + C d/dx(D dT/dx) = E
!!
!!  NOTE: UPON FURTHER REVIEW I THINK THE FORM SOLVED IS ACTUALLY:
!!  dT/dt = A T + B dT/dx + C d/dx(D dT/dx) + E

real(wp), dimension(:), intent(in) :: A,B,C,D,E
real(wp), dimension(:), intent(in) :: Ts
real(wp), intent(in) :: Tsminx1, Tsmaxx1, dt
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
integer, intent(in) ::  flagdirichBottom, flagdirichTop

integer, parameter :: ll=2                   !number of lower diagonals

real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh
integer :: ix1,lx1

real(wp), dimension(size(Ts)) :: TR

real(wp), dimension(size(Ts)) :: TRBDF21D

!> ORGANIZE SIZES AND THERMAL CONDUCTIVITY
lx1=size(Ts)
Dh(1) = 0
Dh(2:lx1) = 0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
!    TR(:)=Ts(:)/dt+E(:)
!! boundaries to be overwritten later...  This is now done for each grid point in a separate statement


!! ## TR HALF STEP:  DEFINE A MATRIX USING BANDED STORAGE

call bottom_bcs(flagdirichBottom, ll, Tsminx1, dx1, M, TR)

!> FIRST INTERIOR GRID POINT
ix1=2

!> ix1-1
M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2 &
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2

!> ix1
M(ll+3,ix1)= 1/(dt/2d0)-A(ix1)/2d0 &
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0

!> ix1+1, super-diag.
M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &
         -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))/2d0
M(ll+1,ix1+2) = 0
TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1) &
  -M(ll+1,ix1+2)*Ts(ix1+2)


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  !! do concurrent OK because only indexing already defined things
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0 &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2d0
  M(ll+3,ix1)=1.0/(dt/2d0)-A(ix1)/2d0 &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0
  M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &    !ix1+1, super-diag.
           -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))/2d0
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
  TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) &
    -M(ll+5,ix1-2)*Ts(ix1-2) &
    -M(ll+4,ix1-1)*Ts(ix1-1) &
    -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0)*Ts(ix1) &
    -M(ll+2,ix1+1)*Ts(ix1+1) &
    -M(ll+1,ix1+2)*Ts(ix1+2)
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2) = 0
M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0 &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))/2d0
M(ll+3,ix1)=1.0/(dt/2d0)-A(ix1)/2d0 &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0
M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0 &        !ix1+1, super-diag.
         -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))/2d0
TR(ix1)=Ts(ix1)/(dt/2d0)+E(ix1) &
  -M(ll+5,ix1-2)*Ts(ix1-2) &
  -M(ll+4,ix1-1)*Ts(ix1-1) &
  -(-A(ix1)/2d0+C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1)/2d0+C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)/2d0)*Ts(ix1) &
  -M(ll+2,ix1+1)*Ts(ix1+1)


call top_bcs(flagdirichTop, ll, lx1, Tsmaxx1, dx1, M, TR)

!! ### TR HALF STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TR,kl=2)



!! ## BDF2 STEP:  DEFINE A MATRIX USING BANDED STORAGE

call bottom_bcs(flagdirichBottom, ll, Tsminx1, dx1, M, TRBDF21D)

!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1.0/(dt/3._wp)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2) = 0
TRBDF21D(ix1)=E(ix1) &
  -1._wp/3._wp*Ts(ix1)/(dt/3._wp) &
  +4._wp/3._wp*TR(ix1)/(dt/3._wp)


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  M(ll+5,ix1-2) = 0                                               !ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &        !ix1-1
             +B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1)=1.0/(dt/3._wp)-A(ix1) &                           !ix1
           +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
           +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &    !ix1+1, super-diag.
           -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2) = 0                                               !ix1+2 grid point
  TRBDF21D(ix1)=E(ix1) &
    -1._wp/3._wp*Ts(ix1)/(dt/3._wp) &
    +4._wp/3._wp*TR(ix1)/(dt/3._wp)
end do


!LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2) = 0
M(ll+4,ix1-1)=-1*C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)=1.0/(dt/3._wp)-A(ix1) &                               !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1)=-1*C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &        !ix1+1, super-diag.
         -1*B(ix1)/(dx1(ix1+1)+dx1(ix1))
TRBDF21D(ix1)=E(ix1) &
  -1._wp/3._wp*Ts(ix1)/(dt/3._wp) &
  +4._wp/3._wp*TR(ix1)/(dt/3._wp)

call top_bcs(flagdirichTop, ll, lx1, Tsmaxx1, dx1, M, TRBDF21D)

!! ## BDF2 STEP MATRIX SOLUTION:  CALL LAPACK'S BANDED SOLVER

!> BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,TRBDF21D,kl=2)

end function TRBDF21D


function backEuler1D(Ts,A,B,C,D,E,Tsminx1,Tsmaxx1,dt,dx1,dx1i, flagdirichBottom, flagdirichTop)

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
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i   !ith centered difference
integer, intent(in) :: flagdirichBottom, flagdirichTop

integer, parameter :: ll=2                   !number of lower diagonals

real(wp), dimension(3*ll+1,size(Ts)) :: M    !note extra rows for lapack workspace
real(wp), dimension(size(Ts)) :: Dh

real(wp), dimension(size(Ts)) :: backEuler1D

integer :: ix1,lx1

!------------------------------------------------------------
!-------DEFINE A MATRIX USING BANDED STORAGE
!------------------------------------------------------------
lx1=size(Ts)
Dh(1) = 0
Dh(2:lx1) = 0.5*(D(1:lx1-1)+D(2:lx1))         !ith left cell wall thermal conductivity
backEuler1D(:)=Ts(:)/dt+E(:)                !boundaries to be overwritten later...

call bottom_bcs(flagdirichBottom, ll, Tsminx1, dx1, M, backEuler1D)

!> FIRST INTERIOR GRID POINT
ix1=2
M(ll+4,ix1-1) = -C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) &            !ix1-1, sub-diaginal
           +B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+3,ix1)= 1/dt-A(ix1) &                                     !ix1
         +C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) &
         +C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
M(ll+2,ix1+1) = -C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) - &        !ix1+1, super-diag.
                B(ix1)/(dx1(ix1+1)+dx1(ix1))
M(ll+1,ix1+2) = 0


!> INTERIOR GRID POINTS
do concurrent (ix1=3:lx1-2)
  M(ll+5,ix1-2) = 0
  !! ix1-2 grid point, sub-diag.
  M(ll+4,ix1-1) = -C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) + &            !ix1-1
                  B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+3,ix1) = 1/dt - A(ix1) + &                                     !ix1
                C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) + &
                C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1)
  M(ll+2,ix1+1) = -C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) - &        !ix1+1, super-diag.
             B(ix1)/(dx1(ix1+1)+dx1(ix1))
  M(ll+1,ix1+2) = 0
  !! ix1+2 grid point
end do


!> LAST INTERIOR GRID POINT
ix1=lx1-1
M(ll+5,ix1-2) = 0
M(ll+4,ix1-1) = -C(ix1)*Dh(ix1)/dx1i(ix1)/dx1(ix1) + &            !ix1-1
                   B(ix1)/(dx1(ix1+1) + dx1(ix1))
M(ll+3,ix1) = 1 / dt - A(ix1) + &                                     !ix1
                C(ix1)*Dh(ix1+1) / dx1i(ix1)/dx1(ix1+1) + &
                C(ix1)*Dh(ix1) / dx1i(ix1) / dx1(ix1)
M(ll+2,ix1+1)= -C(ix1)*Dh(ix1+1)/dx1i(ix1)/dx1(ix1+1) - &        !ix1+1, super-diag.
                 B(ix1) / (dx1(ix1+1) + dx1(ix1))

call top_bcs(flagdirichTop, ll, lx1, Tsmaxx1, dx1, M, backEuler1D)

!! CALL LAPACK'S BANDED SOLVER (INPUT MATRIX MUST BE SHIFTED 'DOWN' BY KL ROWS)
call gbsv(M,backEuler1D,kl=2)

end function backEuler1D


subroutine bottom_bcs(flagdirich, ll, Ts, dx1, M, val)

integer, intent(in) :: flagdirich, ll
real(wp), intent(in) :: Ts, dx1(0:)
real(wp), intent(inout) :: M(:, :), val(:)

integer, parameter :: i = 1

select case (flagdirich)
case (1)
  !! MINX1 BOUNDARY (DIRICHLET)
  M(ll+3,i) = 1
  M(ll+2,i+1) = 0
  M(ll+1,i+2) = 0
  val(i) = Ts
case (0)
  !! if Neumann version, use a 1st order forward...
  M(ll+3,i)=-1.0/dx1(i+1)       !< main diagonal denoted temperature at this grid point... 1*Ts,i=Tsminx1
  M(ll+2,i+1)=1.0/dx1(i+1)      !< 1st super diagonal
  M(ll+1,i+2) = 0                 !< 2nd super diagonal
  val(i) = Ts
  !! here this is not intepreted as temperature, but instead the -1*heat flux divided by thermal conductivity
case default
  write(stderr,*) 'PDEparabolic:bottom_bcs: unknown flagdirich:',flagdirich
  error stop
end select

end subroutine bottom_bcs


subroutine top_bcs(flagdirich, ll, lx1, Ts, dx1, M, val)

integer, intent(in) :: flagdirich, ll, lx1
real(wp), intent(in) :: Ts, dx1(0:)
real(wp), intent(inout) :: M(:, :), val(:)

integer :: i

i = lx1

select case (flagdirich)
case (1)
  !! MAXX1 dirichlet BOUNDARY
  M(ll+5,i-2) = 0
  M(ll+4,i-1) = 0
  M(ll+3,i) = 1
  val(i) = Ts
case (0)
  !!Neumann conditions...
  M(ll+5,i-2) = 0            !< superdiagonal
  M(ll+4,i-1) = -1.0 / dx1(i)  !< subdiagonal
  M(ll+3,i)   = 1.0 / dx1(i)     !< main diag.
  val(i) = Ts     !< here interpreted as -1*heat flux divided by thermal conductivity...
case default
  write(stderr,*) 'PDEparabolic:top_bcs: unknown flagdirich:',flagdirich
  error stop
end select

end subroutine top_bcs

end module PDEparabolic
