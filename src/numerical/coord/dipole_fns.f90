submodule(meshobj_dipole) dipole_fns

!> This submodule contains the functions for which we need to find roots in order to transform
!   dipole to spherical coordinates.

implicit none (type, external)

! options structure for Newton iterations
type(newtopts) :: newtparms

contains

!> convert a single q,p pair into r,theta
module procedure qp2rtheta
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
end procedure qp2rtheta


!> convert a single point r,theta to q,p
module procedure rtheta2qp
  q=Re**2/r**2*cos(theta)
  p=r/Re/(sin(theta)**2)
end procedure rtheta2qp


!> find theta given q,r
module procedure qr2theta
  theta=acos(q*(r/Re)**2)
end procedure qr2theta


!> objective function for newton iterations for solutions of roots for r
module procedure rpoly
  real(wp) ::  q,p

  q=parms(1); p=parms(2);
  fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
end procedure rpoly


!> derivative objective function for newton iterations for roots of r
module procedure rpoly_deriv
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
end procedure rpoly_deriv

end submodule dipole_fns
