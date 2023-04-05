module meshobj_dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems
! Notes:
! 1) may want to overload the calculation procedures (e.g. for metric factors) so that can accept an index or a spherical coordiante input.

! uses
use phys_consts, only: wp,Re,pi,Mmag,mu0,Gconst,Me
use meshobj, only: curvmesh
use newton, only: newtopts,newton_exact,objfun,objfun_deriv
use spherical, only: er_spherical,etheta_spherical,ephi_spherical

implicit none (type, external)
private
public :: dipolemesh, qp2rtheta

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
  real(wp), dimension(:,:,:), pointer :: hqphii,hpphii,hphiphii
  real(wp), dimension(:,:,:,:), pointer :: eq,ep,ephidip
  real(wp), dimension(:,:,:), pointer :: gq,gp,gphi

  contains
    !> Specific methods
    procedure :: calc_rtheta_2D, calc_qp_2D

    !> Bind deferred procedures
    procedure :: init=>init_dipolemesh
    procedure :: make=>make_dipolemesh
    procedure :: calc_er=>calc_er_spher
    procedure :: calc_etheta=>calc_etheta_spher
    procedure :: calc_ephi=>calc_ephi_spher
    procedure :: calc_e1=>calc_eq
    procedure :: calc_e2=>calc_ep
    procedure :: calc_e3=>calc_ephi_dip
    procedure :: calc_grav=>calc_grav_dipole
    procedure :: calc_Bmag=>calc_Bmag_dipole
    procedure :: calc_inclination=>calc_inclination_dipole
    procedure, nopass :: calc_h1=>calc_hq
    procedure, nopass :: calc_h2=>calc_hp
    procedure, nopass :: calc_h3=>calc_hphi_dip
    procedure :: native2ECEFspher=>dipole2ECEFspher

    !> type deallocations, etc.
    final :: destructor
end type dipolemesh


!> declarations and interfaces for submodule functions, apparently these need to be generic interfaces.  These are generally
!   routines that do not directly deal with the derived type data arrays but instead perform very basic calculations
!   related to specifically dipole coordinate transformations.
interface    ! dipole_fns.f90
  module subroutine qp2rtheta(q,p,r,theta)
    real(wp), intent(in) :: q,p
    real(wp), intent(out) :: r,theta
  end subroutine qp2rtheta
  module subroutine rtheta2qp(r,theta,q,p)
    real(wp), intent(in) :: r,theta
    real(wp), intent(out) :: q,p
  end subroutine rtheta2qp
  module elemental function qr2theta(q,r) result(theta)
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
  class(dipolemesh), intent(inout) :: self

  if (.not. self%xi_alloc_status) error stop ' must have curvilinear coordinates defined prior to call init_dipolemesh()'

  ! allocate array space using base type-bound procedure
  call self%calc_coord_diffs()
  call self%init_storage()

  ! now we must associate pointers for extended type alias variables.  This is mostly done in order
  !  to have more readable code.
  self%q=>self%x1; self%p=>self%x2; self%phidip=>self%x3
  self%qint=>self%x1i; self%pint=>self%x2i; self%phiint=>self%x3i
  self%hq=>self%h1; self%hp=>self%h2; self%hphi=>self%h3
  self%hqqi=>self%h1x1i; self%hpqi=>self%h2x1i; self%hphiqi=>self%h3x1i
  self%hqpi=>self%h1x2i; self%hppi=>self%h2x2i; self%hphipi=>self%h3x2i
  self%hqphii=>self%h1x3i; self%hpphii=>self%h2x3i; self%hphiphii=>self%h3x3i
  self%eq=>self%e1; self%ep=>self%e2; self%ephidip=>self%e3
  self%gq=>self%g1; self%gp=>self%g2; self%gphi=>self%g3
end subroutine init_dipolemesh


!> create a dipole mesh structure out of given q,p,phi spacings.  We assume here that the input cell center locations
!   are provide with ghost cells included (note input array indexing in dummy variable declarations.  For new we assume
!   that the fortran code will precompute and store the "full" grid information to save time (but this uses more memory).
subroutine make_dipolemesh(self)
  class(dipolemesh), intent(inout) :: self

  integer :: lqg,lpg,lphig,lq,lp,lphi
  integer :: iphi
  real(wp), dimension(:,:,:), pointer :: r,theta,phispher     ! so these can serve as targets
  real(wp), dimension(:,:,:), pointer :: rqint,thetaqint,phiqint
  real(wp), dimension(:,:,:), pointer :: rpint,thetapint,phipint

  ! check that pointers are correctly associated, which implies that all space has been allocated :)
  if (.not. associated(self%q)) error stop  &
             ' pointers to grid coordiante arrays must be associated prior calling make_dipolemesh()'

  ! size of arrays, including ghost cells
  lqg = size(self%q, 1)
  if(lqg < 1) error stop "lqg must be strictly positive"
  lpg = size(self%p, 1)
  if(lpg < 1) error stop "lpg must be strictly positive"
  lphig = size(self%phidip, 1)
  if(lphig < 1) error stop "lphig must be strictly positive"
  allocate(r(-1:lqg-2, -1:lpg-2, -1:lphig-2), theta(-1:lqg-2, -1:lpg-2, -1:lphig-2))
  allocate(phispher(-1:lqg-2, -1:lpg-2, -1:lphig-2))

! array sizes without ghost cells for convenience
  print "(A,1X,I0,1X,I0,1X,I0)", 'make_dipolemesh:  allocating space for grid of size:',lqg,lpg,lphig
  lq=lqg-4; lp=lpg-4; lphi=lphig-4;
  allocate(rqint(1:lq+1,1:lp,1:lphi))    ! these are just temp vars. needed to compute metric factors
  allocate(thetaqint,phiqint, mold=rqint)
  allocate(rpint(1:lq,1:lp+1,1:lphi))
  allocate(thetapint,phipint, mold=rpint)

  ! convert the cell centers to spherical ECEF coordinates, then tile for longitude dimension
  print*, 'make_dipolemesh:  converting cell centers to spherical coordinates...'
  call self%native2ECEFspher(self%glonctr,self%glatctr,self%q,self%p,self%phidip,r,theta,phispher)

  ! locations of the cell interfaces in q-dimension (along field lines)
  print*, 'make_dipolemesh:  converting cell interfaces in q...'
  call self%native2ECEFspher(self%glonctr,self%glatctr,self%qint,self%p(1:lp),self%phidip(1:lphi),rqint,thetaqint,phiqint)

  ! locations of cell interfaces in p-dimesion (along constant L-shell)
  print*, 'make_dipolemesh:  converting cell interfaces in p...'
  call self%native2ECEFspher(self%glonctr,self%glatctr,self%q(1:lq),self%pint,self%phidip(1:lphi),rpint,thetapint,phipint)

  ! compute and store the metric factors; these need to include ghost cells
  print*, 'make_dipolemesh:  metric factors for cell centers...'
  self%hq(-1:lq+2,-1:lp+2,-1:lphi+2)=self%calc_h1(r,theta,phispher)
  self%hp(-1:lq+2,-1:lp+2,-1:lphi+2)=self%calc_h2(r,theta,phispher)
  self%hphi(-1:lq+2,-1:lp+2,-1:lphi+2)=self%calc_h3(r,theta,phispher)

  ! now assign structure elements and deallocate unneeded temp variables
  !self%r=r(1:lq,1:lp,1:lphi); self%theta=theta(1:lq,1:lp,1:lphi); self%phi=phispher(1:lq,1:lp,1:lphi)   ! don't need ghost cells!
  self%r=r(-1:lq+2,-1:lp+2,-1:lphi+2); self%theta=theta(-1:lq+2,-1:lp+2,-1:lphi+2);    ! keep ghost cells!
  self%phi=phispher(-1:lq+2,-1:lp+2,-1:lphi+2)

  ! compute the geographic coordinates
  print*, 'make_dipolemesh:  geographic coordinates from magnetic...'
  call self%calc_geographic()

  ! q cell interface metric factors
  print*, 'make_dipolemesh:  metric factors for cell q-interfaces...'
  self%hqqi=self%calc_h1(rqint,thetaqint,phispher)
  self%hpqi=self%calc_h2(rqint,thetaqint,phispher)
  self%hphiqi=self%calc_h3(rqint,thetaqint,phispher)

  ! p cell interface metric factors
  print*, 'make_dipolemesh:  metric factors for cell p-intefaces...'
  self%hqpi=self%calc_h1(rpint,thetapint,phispher)
  self%hppi=self%calc_h2(rpint,thetapint,phispher)
  self%hphipi=self%calc_h3(rpint,thetapint,phispher)

  print*, 'make_dipolemesh:  metric factors for cell phi-interfaces...'
  !print*, shape(self%hqphii),shape(self%hpphii),shape(self%hphiphii)
  !print*, shape(self%hq), shape(self%hp), shape(self%hphi)
  self%hqphii(1:lq,1:lp,1:lphi)=self%hq(1:lq,1:lp,1:lphi)         ! note these are not a function of x3 so can just copy things across
  self%hqphii(1:lq,1:lp,lphi+1)=self%hqphii(1:lq,1:lp,lphi)
  self%hpphii(1:lq,1:lp,1:lphi)=self%hp(1:lq,1:lp,1:lphi)         ! seg faults without indices???!!!,  b/c pointers???
  self%hpphii(1:lq,1:lp,lphi+1)=self%hpphii(1:lq,1:lp,lphi)
  self%hphiphii(1:lq,1:lp,1:lphi)=self%hphi(1:lq,1:lp,1:lphi)
  self%hphiphii(1:lq,1:lp,lphi+1)=self%hphiphii(1:lq,1:lp,lphi)

  ! we can now deallocate temp position pointers
  deallocate(r,theta,phispher)
  deallocate(rqint,thetaqint,rpint,thetapint)

  ! spherical ECEF unit vectors (expressed in a Cartesian ECEF basis)
  print*, 'make_dipolemesh:  spherical ECEF unit vectors...'
  call self%calc_er()
  call self%calc_etheta()
  call self%calc_ephi()

  ! dipole coordinate system unit vectors (Cart. ECEF)
  print*, 'make_dipolemesh:  dipole unit vectors...'
  call self%calc_e1()
  call self%calc_e2()
  call self%calc_e3()

  ! magnetic field magnitude
  print*, 'make_dipolemesh:  magnetic fields...'
  call self%calc_Bmag()

  ! gravity components
  print*, 'make_dipolemesh:  gravity...'
  call self%calc_grav()

  ! set the status now that coord. specific calculations are done
  self%coord_alloc_status=.true.

  ! now finish by calling procedures from base type
  print*, 'make_dipolemesh:  base type-bound procedure calls...'
  call self%calc_difflengths()     ! differential lengths (units of m)
  call self%calc_inull()           ! null points (non computational)
  call self%calc_gridflag()        ! compute and store grid type

  ! inclination angle for each field line; awkwardly this must go after gridflag is set...
  print*, 'make_dipolemesh:  inclination angle...'
  call self%calc_inclination()
end subroutine make_dipolemesh


!> convert input dipole coordinates into ECEF spherical coordinates
subroutine dipole2ECEFspher(self,glonctr,glatctr,coord1,coord2,coord3,r,theta,phispher)
  class(dipolemesh), intent(in) :: self
  real(wp) :: glonctr,glatctr     ! unused in this implementation
  real(wp), dimension(:), pointer, intent(in) :: coord1,coord2,coord3
  !real(wp), dimension(:,:,:), intent(inout) :: r,theta,phispher
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord2,1):ubound(coord2,1),lbound(coord3,1):ubound(coord3,1)), &
              intent(inout) :: r,theta,phispher
  integer :: lq,lp,lphi,iphi,iphimin,iphimax
  real(wp), dimension(:), pointer :: q,p,phidip

  q=>coord1; p=>coord2; phidip=>coord3;
  lq=size(q); lp=size(p); lphi=size(phidip);
  iphimin=lbound(phidip,1); iphimax=ubound(phidip,1);

  if (lq /= size(r,1) .or. lp /= size(r,2) .or. lphi /= size(r,3) ) then
    print*, lq,lp,lphi,size(r,1),size(r,2),size(r,3)
    error stop 'dipolemesh::dipole2ECEFspher - r and native coord. arrays not conformable...'
  end if

  call self%calc_rtheta_2D(q,p,r(:,:,iphimin),theta(:,:,iphimin))
  do iphi=iphimin+1,iphimax     ! tile
    r(:,:,iphi)=r(:,:,iphimin)
    theta(:,:,iphi)=theta(:,:,iphimin)
  end do
  do iphi=iphimin,iphimax
    phispher(:,:,iphi)=phidip(iphi)   !scalar assignment should work...
  end do
end subroutine dipole2ECEFspher


subroutine calc_grav_dipole(self)
!! compute gravitational field components
  class(dipolemesh), intent(inout) :: self
  real(wp), dimension(1:self%lx1, 1:self%lx2, 1:self%lx3) :: gr
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  if(any(shape(gr) < 1)) error stop "meshobj_dipole:calc_grav_dipole: lx1,lx2,lx3 must be strictly positive"

  gr = -Gconst*Me / self%r(1:lx1,1:lx2,1:lx3)**2
  !! radial component of gravity
  self%gq = gr*sum(self%er*self%eq, dim=4)
  self%gp = gr*sum(self%er*self%ep, dim=4)
  !gphi=gr*sum(er*ephi,dim=4)
  self%gphi = 0
  !! force to zero to avoid really small values from numerical error
end subroutine calc_grav_dipole


!> compute the magnetic field strength
subroutine calc_Bmag_dipole(self)
  class(dipolemesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%Bmag=mu0*Mmag/4/pi/self%r(1:lx1,1:lx2,1:lx3)**3*sqrt(3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2+1)
end subroutine calc_Bmag_dipole


subroutine calc_inclination_dipole(self)
!! compute the inclination angle (degrees) for each geomagnetic field line
  class(dipolemesh), intent(inout) :: self

  integer :: lq
  real(wp), dimension(1:self%lx1, 1:self%lx2, 1:self%lx3) :: proj

  if(any(shape(proj) < 1)) error stop "meshobj_dipole:calc_inclination_dipole: lx1,lx2,lx3 must be strictly positive"

  lq=size(self%er,1)
  proj=sum(self%er*self%eq,dim=4)
  proj=acos(proj)
  if (self%gridflag==0) then
    !! for a closed grid average over half the domain
    self%I = sum(proj(1:lq/2,:,:), dim=1) / (real(lq, wp)/2)
    !! note use of integer division and casting to real for avging
  else
    !! otherwise average over full domain
    self%I = sum(proj,dim=1) / real(lq, wp)
  end if
  self%I= 90 - min(self%I, pi-self%I) * 180/pi
end subroutine calc_inclination_dipole


!> compute metric factors for q
function calc_hq(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checking
  r=>coord1; theta=>coord2; phi=>coord3;
  hval=r**3/Re**2/sqrt(1+3*cos(theta)**2)
end function calc_hq


!> compute p metric factors
function calc_hp(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checkign

  r=>coord1; theta=>coord2; phi=>coord3;
  hval=Re*sin(theta)**3/sqrt(1+3*cos(theta)**2)
end function calc_hp


!> compute phi metric factor
function calc_hphi_dip(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checking

  r=>coord1; theta=>coord2; phi=>coord3;
  hval=r*sin(theta)
end function calc_hphi_dip


!> radial unit vector (expressed in ECEF cartesian coodinates, components permuted as ix,iy,iz)
subroutine calc_er_spher(self)
  class(dipolemesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%er=er_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_er_spher


!> zenith angle unit vector (expressed in ECEF cartesian coodinates
subroutine calc_etheta_spher(self)
  class(dipolemesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%etheta=etheta_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_etheta_spher


!> azimuth angle unit vector (ECEF cart.)
subroutine calc_ephi_spher(self)
  class(dipolemesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%ephi=ephi_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_ephi_spher


!> unit vector in the q direction
subroutine calc_eq(self)
  class(dipolemesh), intent(inout) :: self
  real(wp), dimension(1:self%lx1,1:self%lx2,1:self%lx3) :: denom
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  denom=sqrt(1+3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2)
  self%eq(:,:,:,1)=-3*cos(self%theta(1:lx1,1:lx2,1:lx3))*sin(self%theta(1:lx1,1:lx2,1:lx3))* &
                        cos(self%phi(1:lx1,1:lx2,1:lx3))/denom
  self%eq(:,:,:,2)=-3*cos(self%theta(1:lx1,1:lx2,1:lx3))*sin(self%theta(1:lx1,1:lx2,1:lx3))* &
                        sin(self%phi(1:lx1,1:lx2,1:lx3))/denom
  self%eq(:,:,:,3)=(1-3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2)/denom   !simplify?
end subroutine calc_eq


!> unit vector in the p direction
subroutine calc_ep(self)
  class(dipolemesh), intent(inout) :: self
  real(wp), dimension(1:self%lx1,1:self%lx2,1:self%lx3) :: denom
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  denom=sqrt(1+3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2)
  self%ep(:,:,:,1)=(1-3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2)*cos(self%phi(1:lx1,1:lx2,1:lx3))/denom
  self%ep(:,:,:,2)=(1-3*cos(self%theta(1:lx1,1:lx2,1:lx3))**2)*sin(self%phi(1:lx1,1:lx2,1:lx3))/denom
  self%ep(:,:,:,3)=3*cos(self%theta(1:lx1,1:lx2,1:lx3))*sin(self%theta(1:lx1,1:lx2,1:lx3))/denom
end subroutine calc_ep


!> unit vector in the phi direction
subroutine calc_ephi_dip(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%ephidip=self%ephi
end subroutine calc_ephi_dip


!> convert a 1D arrays of q,p (assumed to define a 2D grid) into r,theta on a 2D mesh
!   this should be agnostic to the array start index; here just remap as 1:size(array,1), etc.
!   though the dummy argument declarations.  This is necessary due to the way that
subroutine calc_rtheta_2D(self,q,p,r,theta)
  class(dipolemesh), intent(in) :: self
  real(wp), dimension(:), intent(in) :: q
  real(wp), dimension(:), intent(in) :: p
  real(wp), dimension(:,:), intent(inout) :: r,theta

  integer :: iq,ip,lq,lp

  lq = size(self%q)
  !! avoid unused argument warning

  lq=size(q,1); lp=size(p,1);

  do iq=1,lq
    do ip=1,lp
      call qp2rtheta(q(iq),p(ip),r(iq,ip),theta(iq,ip))
    end do
  end do
end subroutine calc_rtheta_2D


!> convert a set of r,theta points (2D arrays) to 2D arrays of q,p
subroutine calc_qp_2D(self,r,theta,q,p)
  class(dipolemesh), intent(in) :: self
  real(wp), dimension(:,:), intent(in) :: r,theta
  real(wp), dimension(:,:), intent(inout) :: q,p

  integer :: i1,i2,ldim1,ldim2

  ldim1 = size(self%r)
  !! avoid unused argument warning

  ldim1=size(r,1); ldim2=size(r,2);

  do i1=1,ldim1
    do i2=1,ldim2
      call rtheta2qp(r(i1,i2),theta(i1,i2),q(i1,i2),p(i1,i2))
    end do
  end do
end subroutine calc_qp_2D


!> type destructor; written generally, viz. as if it is possible some grid pieces are allocated an others are not
subroutine destructor(self)
  type(dipolemesh), intent(inout) :: self

  call self%dissociate_pointers()
  print*, '  dipolemesh destructor completed successfully'
end subroutine destructor

end module meshobj_dipole
