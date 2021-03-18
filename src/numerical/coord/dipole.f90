module dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems

! uses
use newton, only : newton_exact,objfun,objfun_deriv
! use phys_consts, only : wp,Re,pi
integer, parameter :: wp=8
real(wp), parameter :: Re=6370e3
real(wp), parameter :: pi=3.141592


! magnetic pole location
real(wp), private, parameter :: thetan=11*pi/180
real(wp), private, parameter :: phin=289*pi/180

!! procedure declarations
!module procedure(objfun) :: rpoly
!module procedure(objfun_derv) :: rpoly_deriv

contains


subroutine qp2rtheta(q,p,r,theta)
  real(wp), intent(in) :: q,p
  real(wp), intent(out) :: r,theta

end subroutine qp2rtheta


elemental real(8) function qr2theta(q,r) result(theta)
  real(wp), intent(in) :: q,r

  theta=acos(q*(r/Re)**2)
end function qr2theta


subroutine rtheta2qp(r,theta,q,p)
  real(wp), intent(in) :: r,theta
  real(wp), intent(out) :: q,p


end subroutine rtheta2qp


module procedure(objfun) rpoly
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval=q**2*(r/Re)**4 + 1/p*(r/Re) - 1
end procedure rpoly


module procedure(objfun_deriv) rpoly_deriv
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fprimeval=4/Re*q**2*(r/Re)**3 + 1/p/Re
end procedure rpoly_deriv


end module dipole
