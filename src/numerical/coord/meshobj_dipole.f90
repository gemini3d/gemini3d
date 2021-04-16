module meshobj_dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems
! Notes:
! 1) may want to overload the calculation procedures (e.g. for metric factors) so that can accept an index or a spherical coordiante input.

! uses
use phys_consts, only: wp,Re,pi,Mmag,mu0,Gconst,Me
use meshobj, only: curvmesh
use newton, only: newtopts,newton_exact,objfun,objfun_deriv

implicit none (type, external)


! type extension for dipolemesh
type, extends(curvmesh) :: dipolemesh
  real(wp), dimension(:), pointer :: q
  real(wp), dimension(:), pointer :: p
  real(wp), dimension(:), pointer :: phidip
  real(wp), dimension(:), pointer :: qint           ! cell interface locations
  real(wp), dimension(:), pointer :: pint
  real(wp), dimension(:), pointer :: phiint
  real(wp), dimension(:,:,:), pointer :: hq,hp,hphi
  real(wp), dimension(:,:,:), pointer :: hqqi,hpqi,hphiqi
  real(wp), dimension(:,:,:), pointer :: hqpi,hppi,hphipi
  real(wp), dimension(:,:,:,:), pointer :: eq,ep
  real(wp), dimension(:,:,:), pointer :: gq,gp,gphi

  contains

  !> Use generic type-bound procedures here to overload scalar vs. rank 3
  generic :: calc_er=>er_scalar,er_rank3
  procedure :: er_scalar,er_rank3

  generic :: calc_etheta=>etheta_scalar,etheta_rank3
  procedure :: etheta_scalar,etheta_rank3

  generic :: calc_ephi=>ephi_scalar,ephi_rank3
  procedure :: ephi_scalar,ephi_rank3

  generic :: calc_eq=>eq_scalar,eq_rank3
  procedure :: eq_scalar,eq_rank3

  generic :: calc_ep=>ep_scalar,ep_rank3
  procedure :: ep_scalar,ep_rank3

  !> Specific methods
  procedure :: init_dipolemesh
  procedure :: make_dipolemesh
  procedure, private :: calc_grav
  procedure, private :: calc_Bmag
  procedure, private :: calc_inclination
  procedure, private :: calc_hq,calc_hp,calc_hphi
  procedure, private :: calc_geographic
  procedure :: calc_rtheta_2D, calc_qp_2D
end type dipolemesh


!> declarations and interfaces for submodule functions, apparently these need to be generic interfaces.  These are generally
!   routines that do not directly deal with the derived type data arrays but instead perform very basic calculations
!   related to the coordinate transformations.  
interface geomag2geog
  module procedure geomag2geog_scalar
  module procedure geomag2geog_rank3
end interface geomag2geog
interface geog2geomag
  module procedure geog2geomag_scalar
  module procedure geog2geomag_rank3
end interface geog2geomag
interface    ! dipole_fns.f90, spec for submodule functions
  module subroutine geomag2geog_scalar(phi,theta,glon,glat)
    real(wp), intent(in) :: phi,theta
    real(wp), intent(out) :: glon,glat
  end subroutine geomag2geog_scalar
  module subroutine geomag2geog_rank3(phi,theta,glon,glat)
    real(wp), dimension(:,:,:), intent(in) :: phi,theta
    real(wp), dimension(:,:,:), intent(out) :: glon,glat
  end subroutine geomag2geog_rank3
  module subroutine geog2geomag_scalar(glon,glat,phi,theta)
    real(wp), intent(in) :: glon,glat
    real(wp), intent(out) :: phi,theta
  end subroutine geog2geomag_scalar
  module subroutine geog2geomag_rank3(glon,glat,phi,theta)
    real(wp), dimension(:,:,:), intent(in) :: glon,glat
    real(wp), dimension(:,:,:), intent(out) :: phi,theta
  end subroutine geog2geomag_rank3
  elemental module function r2alt(r) result(alt)
    real(wp), intent(in) :: r
    real(wp) :: alt
  end function r2alt
  elemental module function alt2r(alt) result(r)
    real(wp), intent(in) :: alt
    real(wp) :: r
  end function alt2r
  module subroutine qp2rtheta(q,p,r,theta)
    real(wp), intent(in) :: q,p
    real(wp), intent(out) :: r,theta
  end subroutine qp2rtheta
  module subroutine rtheta2qp(r,theta,q,p)
    real(wp), intent(in) :: r,theta
    real(wp), intent(out) :: q,p
  end subroutine rtheta2qp
  elemental module function qr2theta(q,r) result(theta)
    real(wp), intent(in) :: q,r
    real(wp) :: theta
  end function qr2theta
  module function rpoly(x,parms) result(fval)
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
    real(wp) :: fval
  end function rpoly
  module function rpoly_deriv(x,parms) result(fval_deriv)
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
    real(wp) :: fval_deriv
  end function rpoly_deriv
end interface

contains


!> allocate space and associate pointers with arrays in base class.  must runs set_coords first.
subroutine init_dipolemesh(self)
  class(dipolemesh) :: self

  if (.not. self%xi_alloc_status) error stop ' must have curvilinear coordinates defined prior to call init_dipolemesh()'

  ! allocate array space using base type-bound procedure
  call self%calc_coord_diffs()
  call self%init_storage()
  ! fixme: need to add geographic coord arrays first...
  !call self%calc_inull()

  ! now we must associate pointers for extended type alias variables
  self%q=>self%x1; self%p=>self%x2; self%phidip=>self%x3
  self%qint=>self%x1i; self%pint=>self%x2i; self%phiint=>self%x3i
  self%hq=>self%h1; self%hp=>self%h2; self%hphi=>self%h3
  self%hqqi=>self%h1x1i; self%hpqi=>self%h2x1i; self%hphiqi=>self%h3x1i
  self%hqpi=>self%h1x2i; self%hppi=>self%h2x2i; self%hphipi=>self%h3x2i
  ! fixme: add phi interface metric factors
  self%eq=>self%e1; self%ep=>self%e2
  self%gq=>self%g1; self%gp=>self%g2; self%gphi=>self%g3 
end subroutine init_dipolemesh


!> create a dipole mesh structure out of given q,p,phi spacings.  We assume here that the input cell center locations
!   are provide with ghost cells included (note input array indexing in dummy variable declarations.  For new we assume
!   that the fortran code will precompute and store the "full" grid information to save time (but this uses more memory).  
subroutine make_dipolemesh(self) 
  class(dipolemesh) :: self

  integer :: lqg,lpg,lphig,lq,lp,lphi
  integer :: iq,ip,iphi
  real(wp), dimension(:,:,:), allocatable :: r,theta,phispher
  real(wp), dimension(:,:,:), allocatable :: rqint,thetaqint,phiqint
  real(wp), dimension(:,:,:), allocatable :: rpint,thetapint,phipint

  ! check that pointers are correctly associated, which implies that all space has been allocated :)
  if (.not. associated(self%q)) error stop  & 
             ' pointers to grid coordiante arrays must be associated prior calling make_dipolemesh()'

  ! size of arrays, including ghost cells
  lqg=size(self%q,1); lpg=size(self%p,1); lphig=size(self%phidip,1)
  allocate(r(-1:lqg-2,-1:lpg-2,-1:lphig-2),theta(-1:lqg-2,-1:lpg-2,-1:lphig-2))

  ! array sizes without ghost cells for convenience
  print*, ' make_dipolemesh:  allocating space for grid of size:  ',lqg,lpg,lphig
  lq=lqg-4; lp=lpg-4; lphi=lphig-4;
  allocate(phispher(1:lq,1:lp,1:lphi))
  allocate(rqint(1:lq+1,1:lp,1:lphi),thetaqint(1:lq+1,1:lp,1:lphi),phiqint(1:lq+1,1:lp,1:lphi))    ! these are just temp vars. needed to compute metric factors
  allocate(rpint(1:lq+1,1:lp,1:lphi),thetapint(1:lq+1,1:lp,1:lphi),phipint(1:lq+1,1:lp,1:lphi))

  ! convert the cell centers to spherical ECEF coordinates, then tile for longitude dimension
  print*, ' make_dipolemesh:  converting cell centers...'
  call self%calc_rtheta_2D(self%q,self%p,r(:,:,-1),theta(:,:,-1))
  do iphi=0,lphig-2     ! tile
    r(:,:,iphi)=r(:,:,-1)
    theta(:,:,iphi)=theta(:,:,-1)
  end do
  do iphi=1,lphi
    phispher(:,:,iphi)=self%phidip(iphi)   !scalar assignment should work...
  end do

  ! locations of the cell interfaces in q-dimension (along field lines)
  print*, ' make_dipolemesh:  converting cell interfaces in q...'
  call self%calc_rtheta_2D(self%qint,self%p,rqint(:,:,1),thetaqint(:,:,1))
  do iphi=2,lphi
    rqint(:,:,iphi)=rqint(:,:,1)
    thetaqint(:,:,iphi)=thetaqint(:,:,1)
  end do

  ! locations of cell interfaces in p-dimesion (along constant L-shell)
  print*, ' make_dipolemesh:  converting cell interfaces in p...'
  call self%calc_rtheta_2D(self%q,self%pint,rpint(:,:,1),thetapint(:,:,1))
  do iphi=2,lphi
    rpint(:,:,iphi)=rpint(:,:,1)
    thetapint(:,:,iphi)=thetapint(:,:,1)
  end do

  ! now assign structure elements and deallocate unneeded temp variables
  self%r=r(1:lq,1:lp,1:lphi); self%theta=theta(1:lq,1:lp,1:lphi); self%phi=phispher(1:lq,1:lp,1:lphi)
  deallocate(r,theta,phispher)   ! done with these variables now

  ! compute the geographic coordinates
  print*, ' make_dipolemesh:  geographic coordinates from magnetic...'
  call self%calc_geographic(self%r,self%theta,self%phi,self%alt,self%glon,self%glat) 

  ! compute and store the metric factors
  print*, ' make_dipolemesh:  metric factors for cell centers...'
  allocate(self%hq(1:lq,1:lp,1:lphi),self%hp(1:lq,1:lp,1:lphi),self%hphi(1:lq,1:lp,1:lphi))
  self%hq=self%calc_hq(self%r,self%theta)
  self%hp=self%calc_hp(self%r,self%theta)
  self%hphi=self%calc_hphi(self%r,self%theta)

  ! q cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell q-interfaces...'
  self%hqqi=self%calc_hq(rqint,thetaqint)
  self%hpqi=self%calc_hp(rqint,thetaqint)
  self%hphiqi=self%calc_hphi(rqint,thetaqint)

  ! p cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell p intefaces...'
  self%hqpi=self%calc_hq(rpint,thetapint)
  self%hppi=self%calc_hp(rpint,thetapint)
  self%hphipi=self%calc_hphi(rpint,thetapint)

  ! we can now deallocate temp interface arrays
  deallocate(rqint,thetaqint,phiqint,rpint,thetapint,phipint)

  ! spherical ECEF unit vectors (expressed in a Cartesian ECEF basis)
  print*, ' make_dipolemesh:  spherical ECEF unit vectors...'  
  self%er=self%calc_er(self%r,self%theta,self%phi)
  self%etheta=self%calc_etheta(self%r,self%theta,self%phi)
  self%ephi=self%calc_ephi(self%r,self%theta,self%phi)

  ! dipole coordinate system unit vectors (Cart. ECEF)
  print*, ' make_dipolemesh:  dipole unit vectors...'  
  self%eq=self%calc_eq(self%r,self%theta,self%phi)
  self%ep=self%calc_ep(self%r,self%theta,self%phi)
  self%e3=self%ephi

  ! magnetic field magnitude
  print*, ' make_dipolemesh:  magnetic fields...'    
  self%Bmag=self%calc_Bmag(self%r,self%theta)

  ! gravity components
  print*, ' make_dipolemesh:  gravity...'
  call self%calc_grav(self%r,self%eq,self%ep,self%ephi,self%er)

  ! inclination angle for each field line
  print*, ' make_dipolemesh:  inclination angle...'  
  self%I=self%calc_inclination(self%er,self%eq,self%gridflag)

  ! now finish by calling procedures from base type
  self%coord_alloc_status=.true.
  call self%calc_difflengths()     ! differential lengths (units of m)
  call self%calc_inull()
  call self%calc_gridflag()
end subroutine make_dipolemesh


!> compute geographic coordinates of all grid points
subroutine calc_geographic(self,r,theta,phi,alt,glon,glat)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(:,:,:), intent(out) :: alt,glon,glat

  call geomag2geog(phi,theta,glon,glat)
  alt=r2alt(r)
  self%geog_set_status=.true.
end subroutine calc_geographic


!> compute gravitational field components
subroutine calc_grav(self,r,eq,ep,ephi,er)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r
  real(wp), dimension(:,:,:,:), intent(in) :: eq,ep,ephi,er
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: gr
  
  gr=-Gconst*Me/r**2     ! radial component of gravity
  self%gq=gr*sum(self%er*self%eq,dim=4)
  self%gp=gr*sum(self%er*self%ep,dim=4)
  !gphi=gr*sum(er*ephi,dim=4)
  self%gphi=0._wp     ! force to zero to avoid really small values from numerical error
end subroutine calc_grav


!> compute the magnetic field strength
elemental real(wp) function calc_Bmag(self,r,theta) result(Bmag)
  class(dipolemesh), intent(in) :: self
  real(wp), intent(in) :: r,theta

  Bmag=mu0*Mmag/4/pi/r**3*sqrt(3*cos(theta)**2+1)
end function calc_Bmag


!> compute the inclination angle (degrees) for each geomagnetic field line
function calc_inclination(self,er,eq,gridflag) result(Inc)
  class(dipolemesh), intent(in) :: self
  real(wp), dimension(:,:,:,:), intent(in) :: er,eq
  integer, intent(in) :: gridflag
  real(wp), dimension(1:size(er,2),1:size(er,3)) :: Inc
  integer :: lq
  real(wp), dimension(1:size(er,1),1:size(er,2),1:size(er,3)) :: proj

  lq=size(er,1)
  proj=sum(er*eq,dim=4)
  proj=acos(proj)
  if (gridflag==0) then    ! for a closed grid average over half the domain
    Inc=sum(proj(1:lq/2,:,:),dim=1)/(real(lq,wp)/2._wp)    ! note use of integer division and casting to real for avging
  else                     ! otherwise average over full domain
    Inc=sum(proj,dim=1)/real(lq,wp)
  end if
  Inc=90-min(Inc,pi-Inc)*180._wp/pi
end function calc_inclination


!> compute a metric factor for q corresponding to a given r,theta,phi ordered triple
elemental real(wp) function calc_hq(self,r,theta) result(hq)
  class(dipolemesh), intent(in) :: self
  real(wp), intent(in) :: r,theta

  hq=r**3/Re**2/(sqrt(1+3*cos(theta)**2))
end function calc_hq


!> compute p metric factor
elemental real(wp) function calc_hp(self,r,theta) result(hp)
  class(dipolemesh), intent(in) :: self
  real(wp), intent(in) :: r,theta

  hp=Re*sin(theta)**3/(sqrt(1+3*cos(theta)**2))
end function calc_hp


!> compute phi metric factor
elemental real(wp) function calc_hphi(self,r,theta) result(hphi)
  class(dipolemesh), intent(in) :: self
  real(wp), intent(in) :: r,theta

  hphi=r*sin(theta)
end function calc_hphi


!> radial unit vector (expressed in ECEF cartesian coodinates, permuted as ix,iy,iz)
function er_scalar(self,r,theta,phi) result(er)
  class(dipolemesh) :: self
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: er

  er(1)=sin(theta)*cos(phi)
  er(2)=sin(theta)*sin(phi)
  er(3)=cos(theta)
end function er_scalar
function er_rank3(self,r,theta,phi) result(er)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: er

  er(:,:,:,1)=sin(theta)*cos(phi)
  er(:,:,:,2)=sin(theta)*sin(phi)
  er(:,:,:,3)=cos(theta)
end function er_rank3


!> zenith angle unit vector (expressed in ECEF cartesian coodinates
function etheta_scalar(self,r,theta,phi) result(etheta)
  class(dipolemesh) :: self
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: etheta

  etheta(1)=cos(theta)*cos(phi)
  etheta(2)=cos(theta)*sin(phi)
  etheta(3)=-sin(theta)
end function etheta_scalar
function etheta_rank3(self,r,theta,phi) result(etheta)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: etheta

  etheta(:,:,:,1)=cos(theta)*cos(phi)
  etheta(:,:,:,2)=cos(theta)*sin(phi)
  etheta(:,:,:,3)=-sin(theta)
end function etheta_rank3


!> azimuth angle unit vector (ECEF cart.)
function ephi_scalar(self,r,theta,phi) result(ephi)
  class(dipolemesh) :: self
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: ephi

  ephi(1)=-sin(phi)
  ephi(2)=cos(phi)
  ephi(3)=0
end function ephi_scalar
function ephi_rank3(self,r,theta,phi) result(ephi)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: ephi

  ephi(:,:,:,1)=-sin(phi)
  ephi(:,:,:,2)=cos(phi)
  ephi(:,:,:,3)=0
end function ephi_rank3


!> unit vector in the q direction
function eq_scalar(self,r,theta,phi) result(eq)
  class(dipolemesh) :: self
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: eq
  real(wp) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  eq(1)=-3*cos(theta)*sin(theta)*cos(phi)/denom
  eq(2)=-3*cos(theta)*sin(theta)*sin(phi)/denom
  eq(3)=(1-3*cos(theta)**2)/denom   !simplify?
end function eq_scalar
function eq_rank3(self,r,theta,phi) result(eq)
  class(dipolemesh) :: self
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3),3) :: eq
  real(wp), dimension(1:size(r,1),1:size(r,2),1:size(r,3)) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  eq(:,:,:,1)=-3*cos(theta)*sin(theta)*cos(phi)/denom
  eq(:,:,:,2)=-3*cos(theta)*sin(theta)*sin(phi)/denom
  eq(:,:,:,3)=(1-3*cos(theta)**2)/denom   !simplify?
end function eq_rank3


!> unit vector in the p direction
function ep_scalar(self,r,theta,phi) result(ep)
  class(dipolemesh) :: self
  real(wp), intent(in) :: r,theta,phi
  real(wp), dimension(3) :: ep
  real(wp) :: denom  

  denom=sqrt(1+3*cos(theta)**2)
  ep(1)=(1-3*cos(theta)**2)*cos(phi)/denom
  ep(2)=(1-3*cos(theta)**2)*sin(phi)/denom
  ep(3)=3*cos(theta)*sin(theta)/denom
end function ep_scalar
function ep_rank3(self,r,theta,phi) result(ep)
  class(dipolemesh) :: self
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
subroutine calc_rtheta_2D(self,q,p,r,theta)
  class(dipolemesh) :: self
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
end subroutine calc_rtheta_2D


!> convert a set of r,theta points (2D arrays) to 2D arrays of q,p
subroutine calc_qp_2D(self,r,theta,q,p)
  class(dipolemesh) :: self
  real(wp), dimension(:,:), intent(in) :: r,theta
  real(wp), dimension(:,:), intent(out) :: q,p

  integer :: i1,i2,ldim1,ldim2

  ldim1=size(r,1); ldim2=size(r,2);  

  do i1=1,ldim1
    do i2=1,ldim2
      call rtheta2qp(r(i1,i2),theta(i1,i2),q(i1,i2),p(i1,i2))
    end do
  end do
end subroutine calc_qp_2D

end module meshobj_dipole

