module dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems
! Notes:
! 1) may want to overload the calculation procedures (e.g. for metric factors) so that can accept an index or a spherical coordiante input.

! uses
use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use newton, only : newtopts,newton_exact,objfun,objfun_deriv

implicit none (type, external)

! use phys_consts, only : wp,Re,pi,Mmag
real(wp), parameter ::  Re = 6371.0e3_wp
real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
real(wp), parameter :: Mmag=7.94e22
real(wp), parameter :: mu0=4*pi*1e-7
real(wp), parameter :: Gconst=6.67408e-11_wp     !Universal gravitation constant
real(wp), parameter :: Me = 5.9722e24_wp

! magnetic pole location
real(wp), private, parameter :: thetan=11*pi/180
real(wp), private, parameter :: phin=289*pi/180

type(newtopts), public :: newtparms

!> these to be dealt with later, once things are done
!public :: get_rtheta_2D,get_qp_2D,qp2rtheta,rtheta2qp,rpoly,rpoly_deriv,hq
!private :: qr2theta

!> overload interfaces for unit vectors since we can't handle these with elementals...
!   The intent here is to be able to call these for a single point OR the entire grid.
interface er
  module procedure er_scalar
  module procedure er_rank3
end interface
interface etheta
  module procedure etheta_scalar
  module procedure etheta_rank3
end interface etheta
interface ephi
  module procedure ephi_scalar
  module procedure ephi_rank3
end interface ephi
interface eq
  module procedure eq_scalar
  module procedure eq_rank3
end interface eq
interface ep
  module procedure ep_scalar
  module procedure ep_rank3
end interface ep

contains

!> compute gravitational field components
subroutine grav(r,eq,ep,ephi,er,gq,gp,gphi)
  real(wp), dimension(:,:,:), intent(in) :: r
  real(wp), dimension(:,:,:,:), intent(in) :: eq,ep,ephi,er
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)), intent(out) :: gq,gp,gphi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: gr
  
  gr=-Gconst*Me/r**2     ! radial component of gravity
  gq=gr*sum(er*eq,dim=4)
  gp=gr*sum(er*ep,dim=4)
  !gphi=gr*sum(er*ephi,dim=4)
  gphi=0._wp     ! force to zero to avoid really small values from numerical error
end subroutine grav


!> compute the magnetic field strength
elemental real(wp) function Bmag(r,theta)
  real(wp), intent(in) :: r,theta

  Bmag=mu0*Mmag/4/pi/r**3*sqrt(3*cos(theta)**2+1)
end function Bmag

!> compute the inclination angle for each geomagnetic field line


!> compute a metric factor for q corresponding to a given r,theta,phi ordered triple
elemental real(wp) function hq(r,theta,phi)
  real(wp), intent(in) :: r,theta,phi

  hq=r**3/Re**2/(sqrt(1+3*cos(theta)**2))
end function hq

!> compute p metric factor
elemental real(wp) function hp(r,theta,phi)
  real(wp), intent(in) :: r,theta,phi

  hp=Re*sin(theta)**3/(sqrt(1+3*cos(theta)**2))
end function hp

!> compute phi metric factor
elemental real(wp) function hphi(r,theta,phi)
  real(wp), intent(in) :: r,theta,phi

  hphi=r*sin(theta)
end function hphi


!> radial unit vector (expressed in ECEF cartesian coodinates, permuted as ix,iy,iz)
function er_scalar(r,theta,phi) result(er)
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: er

  er(1)=sin(theta)*cos(phi)
  er(2)=sin(theta)*sin(phi)
  er(3)=cos(theta)
end function er_scalar
function er_rank3(r,theta,phi) result(er)
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: er

  er(:,:,:,1)=sin(theta)*cos(phi)
  er(:,:,:,2)=sin(theta)*sin(phi)
  er(:,:,:,3)=cos(theta)
end function er_rank3


!> zenith angle unit vector (expressed in ECEF cartesian coodinates
function etheta_scalar(r,theta,phi) result(etheta)
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: etheta

  etheta(1)=cos(theta)*cos(phi)
  etheta(2)=cos(theta)*sin(phi)
  etheta(3)=-sin(theta)
end function etheta_scalar
function etheta_rank3(r,theta,phi) result(etheta)
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: etheta

  etheta(:,:,:,1)=cos(theta)*cos(phi)
  etheta(:,:,:,2)=cos(theta)*sin(phi)
  etheta(:,:,:,3)=-sin(theta)
end function etheta_rank3


!> azimuth angle unit vector (ECEF cart.)
function ephi_scalar(r,theta,phi) result(ephi)
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: ephi

  ephi(1)=-sin(phi)
  ephi(2)=cos(phi)
  ephi(3)=0
end function ephi_scalar
function ephi_rank3(r,theta,phi) result(ephi)
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: ephi

  ephi(:,:,:,1)=-sin(phi)
  ephi(:,:,:,2)=cos(phi)
  ephi(:,:,:,3)=0
end function ephi_rank3


!> unit vector in the q direction
function eq_scalar(r,theta,phi) result(eq)
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: eq
  real(wp) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  eq(1)=-3*cos(theta)*sin(theta)*cos(phi)/denom
  eq(2)=-3*cos(theta)*sin(theta)*sin(phi)/denom
  eq(3)=(1-3*cos(theta)**2)/denom   !simplify?
end function eq_scalar
function eq_rank3(r,theta,phi) result(eq)
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: eq
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  eq(:,:,:,1)=-3*cos(theta)*sin(theta)*cos(phi)/denom
  eq(:,:,:,2)=-3*cos(theta)*sin(theta)*sin(phi)/denom
  eq(:,:,:,3)=(1-3*cos(theta)**2)/denom   !simplify?
end function eq_rank3


!> unit vector in the p direction
function ep_scalar(r,theta,phi) result(ep)
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: ep
  real(wp) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  ep(1)=(1-3*cos(theta)**2)*cos(phi)/denom
  ep(2)=(1-3*cos(theta)**2)*sin(phi)/denom
  ep(3)=3*cos(theta)*sin(theta)/denom
end function ep_scalar
function ep_rank3(r,theta,phi) result(ep)
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: ep
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  ep(:,:,:,1)=(1-3*cos(theta)**2)*cos(phi)/denom
  ep(:,:,:,2)=(1-3*cos(theta)**2)*sin(phi)/denom
  ep(:,:,:,3)=3*cos(theta)*sin(theta)/denom
end function ep_rank3


!> convert a 1D arrays of q,p (assumed to define a 2D grid) into r,theta on a 2D mesh
subroutine get_rtheta_2D(q,p,r,theta)
  real(wp), dimension(:), intent(in) :: q
  real(wp), dimension(:), intent(in) :: p
  real(wp), dimension(:,:), intent(out) :: r,theta

  integer :: iq,ip,lq,lp

  lq=size(q,1); lp=size(p,1);

  do iq=1,lq
    do ip=1,lp
      call qp2rtheta(q(iq),p(ip),r(iq,ip),theta(iq,ip))
    end do
  end do
end subroutine get_rtheta_2D


!> convert a set of r,theta points (2D arrays) to 2D arrays of q,p
subroutine get_qp_2D(r,theta,q,p)
  real(wp), dimension(:,:), intent(in) :: r,theta
  real(wp), dimension(:,:), intent(out) :: q,p

  integer :: i1,i2,ldim1,ldim2

  ldim1=size(r,1); ldim2=size(r,2);  

  do i1=1,ldim1
    do i2=1,ldim2
      call rtheta2qp(r(i1,i2),theta(i1,i2),q(i1,i2),p(i1,i2))
    end do
  end do
end subroutine get_qp_2D


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
elemental real(wp) function qr2theta(q,r) result(theta)
  real(wp), intent(in) :: q,r

  theta=acos(q*(r/Re)**2)
end function qr2theta


!> objective function for newton iterations for solutions of roots for r
real(wp) function rpoly(x,parms) result(fval)
  real(wp), intent(in) :: x
  real(wp), dimension(:), intent(in) :: parms
  real(wp) ::  q,p

  q=parms(1); p=parms(2);
  fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
end function rpoly


!> derivative objective function for newton iterations for roots of r
real(wp) function rpoly_deriv(x,parms) result(fval_deriv)
  real(wp), intent(in) :: x
  real(wp), dimension(:), intent(in) :: parms
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
end function rpoly_deriv

end module dipole
