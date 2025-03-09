module dipole

!> This submodule contains the functions for which we need to find roots in order to transform
!   dipole to spherical coordinates.

use phys_consts, only: wp, Re
use newton, only: newtopts, objfun, objfun_deriv, newton_exact

implicit none (type, external)

! options structure for Newton iterations
type(newtopts) :: newtparms

private
public :: qp2rtheta, rtheta2qp

contains
  !> convert a single q,p pair into r,theta
  subroutine qp2rtheta(q,p,r,theta)
    real(wp), intent(in) :: q,p
    real(wp), intent(out) :: r,theta
  
    real(wp), dimension(2) :: parms
    real(wp) :: r0
    procedure(objfun), pointer :: f
    procedure(objfun_deriv), pointer :: fprime
    integer :: maxrestart, maxr, r0step
    integer :: it,ir0
    logical :: converged
  
    ! Set parameters of the restart and Newton iterations
    maxrestart=400
    maxr=100*Re
    r0step=0.25*Re
    newtparms%maxit=100
    newtparms%derivtol=1e-18
    newtparms%tol=1e-11
    newtparms%verbose=.false.
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
  
    ! Once we have r can algebraically solve for theta
    theta=qr2theta(q,r)
  end subroutine qp2rtheta
  
  
  !> convert a single point r,theta to q,p
  subroutine rtheta2qp(r,theta,q,p)
    real(wp), intent(in) :: r,theta
    real(wp), intent(out) :: q,p
  
    q=Re**2/r**2*cos(theta)
    p=r/Re/(sin(theta)**2)
  end subroutine rtheta2qp
  
  
  !> find theta given q,r
  elemental function qr2theta(q,r) result(theta)
    real(wp), intent(in) :: q,r
    real(wp) :: theta
  
    theta=acos(q*(r/Re)**2)
  end function qr2theta
  
  
  !> objective function for newton iterations for solutions of roots for r
  function rpoly(x,parms) result(fval)
      real(wp), intent(in) :: x
      real(wp), dimension(:), intent(in) :: parms
      real(wp) :: fval
  
    real(wp) ::  q,p
  
    q=parms(1); p=parms(2);
    fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
  end function rpoly
  
  
  !> derivative objective function for newton iterations for roots of r
  function rpoly_deriv(x,parms) result(fval_deriv)
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
    real(wp) :: fval_deriv
  
    real(wp) :: q,p
  
    q=parms(1); p=parms(2);
    fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
  end function rpoly_deriv
end module dipole
