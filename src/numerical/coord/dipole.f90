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

! options structure for Newton iterations
type(newtopts), public :: newtparms

! class definition for dipolemesh
type :: dipolemesh
  integer :: gridflag

  real(wp), dimension(:), allocatable :: q
  real(wp), dimension(:), allocatable :: p
  real(wp), dimension(:), allocatable :: phi
  real(wp), dimension(:,:,:), allocatable :: r,theta,phispher    ! spherical ECEF coordinates for each grid point
  real(wp), dimension(:), allocatable :: qint           ! cell interface locations
  real(wp), dimension(:), allocatable :: pint
  real(wp), dimension(:), allocatable :: phiint
  real(wp), dimension(:,:,:), allocatable :: hq,hp,hphi
  real(wp), dimension(:,:,:), allocatable :: hqqi,hpqi,hphiqi
  real(wp), dimension(:,:,:), allocatable :: hqpi,hppi,hphipi
  real(wp), dimension(:,:,:,:), allocatable :: er,etheta,ephi
  real(wp), dimension(:,:,:,:), allocatable :: eq,ep
  real(wp), dimension(:,:,:), allocatable :: Bmag
  real(wp), dimension(:,:,:), allocatable :: gq,gp,gphi
  real(wp), dimension(:,:), allocatable :: Inc
end type dipolemesh

!> these to be dealt with later, once things are done
!public :: get_rtheta_2D,get_qp_2D,qp2rtheta,rtheta2qp,rpoly,rpoly_deriv,hq
!private :: qr2theta

!> overload interfaces for unit vectors since we can't handle these with elementals...
!   The intent here is to be able to call these for a single point OR the entire grid.
interface get_er
  module procedure er_scalar
  module procedure er_rank3
end interface get_er
interface get_etheta
  module procedure etheta_scalar
  module procedure etheta_rank3
end interface get_etheta
interface get_ephi
  module procedure ephi_scalar
  module procedure ephi_rank3
end interface get_ephi
interface get_eq
  module procedure eq_scalar
  module procedure eq_rank3
end interface get_eq
interface get_ep
  module procedure ep_scalar
  module procedure ep_rank3
end interface get_ep

contains


!> create a dipole mesh structure out of given q,p,phi spacings.  We assume here that the input cell center locations
!   are provide with ghost cells included (note input array indexing in dummy variable declarations.  For new we assume
!   that the fortran code will precompute and store the "full" grid information to save time (but this uses more memory).  
function make_dipolemesh(q,p,phi) result(x)
  real(wp), dimension(-1:), intent(in) :: q
  real(wp), dimension(-1:), intent(in) :: p
  real(wp), dimension(-1:), intent(in) :: phi
  type(dipolemesh) :: x

  integer :: lqg,lpg,lphig,lq,lp,lphi
  integer :: iq,ip,iphi
  real(wp), dimension(:,:,:), allocatable :: r,theta,phispher
  real(wp), dimension(:,:,:), allocatable :: rqint,thetaqint,phiqint
  real(wp), dimension(:,:,:), allocatable :: rpint,thetapint,phipint

  ! size of arrays, including ghost cells
  lqg=size(q,1); lpg=size(p,1); lphig=size(phi,1)
  allocate(x%q(-1:lqg-2),x%p(-1:lpg-2),x%phi(-1:lphig-2))
  x%q=q; x%p=p; x%phi=phi;
  allocate(r(-1:lqg-2,-1:lpg-2,-1:lphig-2),theta(-1:lqg-2,-1:lpg-2,-1:lphig-2))

  ! array sizes without ghost cells for convenience
  print*, ' make_dipolemesh:  allocating space for grid of size:  ',lqg,lpg,lphig
  lq=lqg-4; lp=lpg-4; lphi=lphig-4;
  allocate(x%qint(1:lq+1),x%pint(1:lp+1),x%phiint(1:lphi+1))     !qint(iq) references the left cell wall location of the ith grid point
  allocate(phispher(1:lq,1:lp,1:lphi))
  allocate(rqint(1:lq+1,1:lp,1:lphi),thetaqint(1:lq+1,1:lp,1:lphi),phiqint(1:lq+1,1:lp,1:lphi))    ! these are just temp vars. needed to compute metric factors
  allocate(rpint(1:lq+1,1:lp,1:lphi),thetapint(1:lq+1,1:lp,1:lphi),phipint(1:lq+1,1:lp,1:lphi))

  ! convert the cell centers to spherical ECEF coordinates, then tile for longitude dimension
  print*, ' make_dipolemesh:  converting cell centers...'
  call get_rtheta_2D(q,p,r(:,:,-1),theta(:,:,-1))
  do iphi=0,lphig-2     ! tile
    r(:,:,iphi)=r(:,:,-1)
    theta(:,:,iphi)=theta(:,:,-1)
  end do
  do iphi=1,lphi
    phispher(:,:,iphi)=phi(iphi)   !scalar assignment should work...
  end do

  ! compute interface coordinates, both dipole and spherical, can be made directly in the mesh
  !  objects memory space
  x%qint(1:lq+1)=1._wp/2._wp*(q(0:lq)+q(1:lq+1))
  x%pint(1:lq+1)=1._wp/2._wp*(p(0:lp)+p(1:lp+1))
  x%phiint(1:lq+1)=1._wp/2._wp*(phi(0:lphi)+phi(1:lphi+1))

  ! locations of the cell interfaces in q-dimension (along field lines)
  print*, ' make_dipolemesh:  converting cell interfaces in q...'
  call get_rtheta_2D(x%qint,p,rqint(:,:,1),thetaqint(:,:,1))
  do iphi=2,lphi
    rqint(:,:,iphi)=rqint(:,:,1)
    thetaqint(:,:,iphi)=thetaqint(:,:,1)
  end do

  ! locations of cell interfaces in p-dimesion (along constant L-shell)
  print*, ' make_dipolemesh:  converting cell interfaces in p...'
  call get_rtheta_2D(q,x%pint,rpint(:,:,1),thetapint(:,:,1))
  do iphi=2,lphi
    rpint(:,:,iphi)=rpint(:,:,1)
    thetapint(:,:,iphi)=thetapint(:,:,1)
  end do

  ! now assign structure elements and deallocate unneeded temp variables
  allocate(x%r(1:lq,1:lp,1:lphi),x%theta(1:lq,1:lp,1:lphi),x%phispher(1:lq,1:lp,1:lphi))
  x%r=r(1:lq,1:lp,1:lphi); x%theta=theta(1:lq,1:lp,1:lphi); x%phispher=phispher(1:lq,1:lp,1:lphi)
  deallocate(r,theta,phispher)

  ! compute and store the metric factors
  print*, ' make_dipolemesh:  metric factors for cell centers...'
  allocate(x%hq(1:lq,1:lp,1:lphi),x%hp(1:lq,1:lp,1:lphi),x%hphi(1:lq,1:lp,1:lphi))
  x%hq=get_hq(x%r,x%theta)
  x%hp=get_hp(x%r,x%theta)
  x%hphi=get_hphi(x%r,x%theta)

  ! q cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell q-interfaces...'
  allocate(x%hqqi(1:lq+1,1:lp,1:lphi),x%hpqi(1:lq+1,1:lp,1:lphi),x%hphiqi(1:lq+1,1:lp,1:lphi))
  x%hqqi=get_hq(rqint,thetaqint)
  x%hpqi=get_hp(rqint,thetaqint)
  x%hphiqi=get_hphi(rqint,thetaqint)

  ! p cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell p intefaces...'
  allocate(x%hqpi(1:lq,1:lp,1+1:lphi),x%hppi(1:lq,1:lp+1,1:lphi),x%hphipi(1:lq,1:lp+1,1:lphi))  
  x%hqpi=get_hq(rpint,thetapint)
  x%hppi=get_hp(rpint,thetapint)
  x%hphipi=get_hphi(rpint,thetapint)

  ! we can now deallocate temp interface arrays
  deallocate(rqint,thetaqint,phiqint,rpint,thetapint,phipint)

  ! spherical ECEF unit vectors (expressed in a Cartesian ECEF basis)
  print*, ' make_dipolemesh:  spherical ECEF unit vectors...'  
  allocate(x%er(1:lq,1:lp,1:lphi,3),x%etheta(1:lq,1:lp,1:lphi,3),x%ephi(1:lq,1:lp,1:lphi,3))
  x%er=get_er(x%r,x%theta,x%phispher)
  x%etheta=get_etheta(x%r,x%theta,x%phispher)
  x%ephi=get_ephi(x%r,x%theta,x%phispher)

  ! dipole coordinate system unit vectors (Cart. ECEF)
  print*, ' make_dipolemesh:  dipole unit vectors...'  
  allocate(x%eq(1:lq,1:lp,1:lphi,3),x%ep(1:lq,1:lp,1:lphi,3))
  x%eq=get_eq(x%r,x%theta,x%phispher)
  x%ep=get_ep(x%r,x%theta,x%phispher)

  ! magnetic field magnitude
  print*, ' make_dipolemesh:  magnetic fields...'    
  allocate(x%Bmag(1:lq,1:lp,1:lphi))
  x%Bmag=get_Bmag(x%r,x%theta)

  ! gravity components
  print*, ' make_dipolemesh:  gravity...'
  allocate(x%gq(1:lq,1:lp,1:lphi),x%gp(1:lq,1:lp,1:lphi),x%gphi(1:lq,1:lp,1:lphi))
  call get_grav(x%r,x%eq,x%ep,x%ephi,x%er,x%gq,x%gp,x%gphi)

  ! FIXME: hardcode grid type for now
  x%gridflag=0

  ! inclination angle for each field line
  print*, ' make_dipolemesh:  inclination angle...'  
  allocate(x%Inc(1:lp,1:lphi))
  x%Inc=get_inclination(x%er,x%eq,x%gridflag)
end function make_dipolemesh


!> deallocate space used by mesh (to be called at the end of the program, presumably)
subroutine de_dipolemesh(x)
  type(dipolemesh), intent(inout) :: x

  deallocate(x%q,x%p,x%phi)
  deallocate(x%r,x%theta,x%phispher)
  deallocate(x%qint,x%pint,x%phiint)
  deallocate(x%hq,x%hp,x%hphi)
  deallocate(x%hqqi,x%hpqi,x%hphiqi)
  deallocate(x%hqpi,x%hppi,x%hphipi)
  deallocate(x%er,x%etheta,x%ephi)
  deallocate(x%eq,x%ep)
  deallocate(x%Bmag)
  deallocate(x%gq,x%gp,x%gphi)
  deallocate(x%Inc)

end subroutine de_dipolemesh


!> compute gravitational field components
subroutine get_grav(r,eq,ep,ephi,er,gq,gp,gphi)
  real(wp), dimension(:,:,:), intent(in) :: r
  real(wp), dimension(:,:,:,:), intent(in) :: eq,ep,ephi,er
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)), intent(out) :: gq,gp,gphi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: gr
  
  gr=-Gconst*Me/r**2     ! radial component of gravity
  gq=gr*sum(er*eq,dim=4)
  gp=gr*sum(er*ep,dim=4)
  !gphi=gr*sum(er*ephi,dim=4)
  gphi=0._wp     ! force to zero to avoid really small values from numerical error
end subroutine get_grav


!> compute the magnetic field strength
elemental real(wp) function get_Bmag(r,theta) result(Bmag)
  real(wp), intent(in) :: r,theta

  Bmag=mu0*Mmag/4/pi/r**3*sqrt(3*cos(theta)**2+1)
end function get_Bmag


!> compute the inclination angle (degrees) for each geomagnetic field line
function get_inclination(er,eq,gridflag) result(Inc)
  real(wp), dimension(:,:,:,:), intent(in) :: er,eq
  integer, intent(in) :: gridflag
  real(wp), dimension(1:size(er,2),1:size(er,3)) :: Inc
  integer :: lq
  real(wp), dimension(1:size(er,1),1:size(er,2),1:size(er,3)) :: proj

  lq=size(er,1)
  proj=sum(er*eq,dim=4)
  if (gridflag==0) then    ! for a closed grid average over half the domain
    Inc=sum(proj(1:lq/2,:,:),dim=1)/real(lq/2,wp)    ! note use of integer division and casting to real for avging
  else                     ! otherwise average over full domain
    Inc=sum(proj,dim=1)/real(lq,wp)
  end if
  Inc=90-min(Inc,pi-Inc)*180._wp/pi
end function get_inclination


!> compute a metric factor for q corresponding to a given r,theta,phi ordered triple
elemental real(wp) function get_hq(r,theta) result(hq)
  real(wp), intent(in) :: r,theta

  hq=r**3/Re**2/(sqrt(1+3*cos(theta)**2))
end function get_hq


!> compute p metric factor
elemental real(wp) function get_hp(r,theta) result(hp)
  real(wp), intent(in) :: r,theta

  hp=Re*sin(theta)**3/(sqrt(1+3*cos(theta)**2))
end function get_hp


!> compute phi metric factor
elemental real(wp) function get_hphi(r,theta) result(hphi)
  real(wp), intent(in) :: r,theta

  hphi=r*sin(theta)
end function get_hphi


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
!   this should be agnostic to the array start index; here just remap as 1:size(array,1), etc.
!   though the dummy argument declarations.  This is necessary due to the way that 
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
