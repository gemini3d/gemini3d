module dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems

! uses
use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use newton, only : newtopts,newton_exact,objfun,objfun_deriv

implicit none (type, external)

! use phys_consts, only : wp,Re,pi
real(wp), parameter :: Re=6370e3
real(wp), parameter :: pi=3.141592

! magnetic pole location
real(wp), private, parameter :: thetan=11*pi/180
real(wp), private, parameter :: phin=289*pi/180

! print output
logical :: verbose=.true.

contains

subroutine qp2rtheta(q,p,r,theta)
  real(wp), intent(in) :: q,p
  real(wp), intent(out) :: r,theta

  real(wp), dimension(2) :: parms
  real(wp) :: r0
  procedure(objfun), pointer :: f
  procedure(objfun_deriv), pointer :: fprime
  integer :: maxrestart, maxr, r0step
  type(newtopts) :: newtparms
  integer :: it,ir0
  logical :: converged


  ! Set parameters of the restart and Newton iterations
  maxrestart=400
  maxr=100*Re
  r0step=0.25*Re
  newtparms%maxit=100
  newtparms%derivtol=1e-18
  newtparms%tol=1e-11
  newtparms%verbose=.true.
  f=>rpoly
  fprime=>rpoly_deriv
  parms=[q,p]

  ! Newton iterations with restarting (see parameters above for limits) until we get a satisfactory result
  r=0; converged=.false.; ir0=1;
  do while (.not. converged .and. ir0<maxrestart .and. (r<=0 .or. r>maxr))
    r0=(ir0-1)*(r0step)    ! change starting point in increments of 0.25 Re until we get a "good" answer
    call newton_exact(f,fprime,r0,parms,newtparms,r,it,converged)
    ir0=ir0+1
  end do
  theta=qr2theta(q,r)

end subroutine qp2rtheta


elemental real(wp) function qr2theta(q,r) result(theta)
  real(wp), intent(in) :: q,r

  theta=acos(q*(r/Re)**2)
end function qr2theta


real(wp) function rpoly(x,parms) result(fval)
  real(wp), intent(in) :: x
  real(wp), dimension(:), intent(in) :: parms
  real(wp) ::  q,p

  q=parms(1); p=parms(2);
  fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
end function rpoly


real(wp) function rpoly_deriv(x,parms) result(fval_deriv)
  real(wp), intent(in) :: x
  real(wp), dimension(:), intent(in) :: parms
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
end function rpoly_deriv

end module dipole
