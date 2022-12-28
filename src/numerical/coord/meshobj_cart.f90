module meshobj_cart

!> Contains data and subroutines for managing a Cartesian mesh

! uses
use phys_consts, only: wp,Re,pi,Mmag,mu0,Gconst,Me
use meshobj, only: curvmesh
use spherical, only: er_spherical,etheta_spherical,ephi_spherical
use geomagnetic, only: geog2geomag,geomag2geog,r2alt,alt2r

implicit none (type, external)
private
public :: cartmesh

! type extension for cartmesh
type, extends(curvmesh) :: cartmesh
  real(wp), dimension(:), pointer :: z
  real(wp), dimension(:), pointer :: x
  real(wp), dimension(:), pointer :: y
  real(wp), dimension(:), pointer :: zint           ! cell interface locations
  real(wp), dimension(:), pointer :: xint
  real(wp), dimension(:), pointer :: yint
  real(wp), dimension(:,:,:), pointer :: hz,hx,hy
  real(wp), dimension(:,:,:), pointer :: hzzi,hxzi,hyzi
  real(wp), dimension(:,:,:), pointer :: hzxi,hxxi,hyxi
  real(wp), dimension(:,:,:), pointer :: hzyi,hxyi,hyyi
  real(wp), dimension(:,:,:,:), pointer :: ez,ex,ey
  real(wp), dimension(:,:,:), pointer :: gz,gx,gy

  contains
    !> Bind deferred procedures
    procedure :: init=>init_cartmesh
    procedure :: make=>make_cartmesh
    procedure :: calc_er=>calc_er_spher
    procedure :: calc_etheta=>calc_etheta_spher
    procedure :: calc_ephi=>calc_ephi_spher
    procedure :: calc_e1=>calc_ez
    procedure :: calc_e2=>calc_ex
    procedure :: calc_e3=>calc_ey
    procedure :: calc_grav=>calc_grav_cart
    procedure :: calc_Bmag=>calc_Bmag_cart
    procedure :: calc_inclination=>calc_inclination_cart
    procedure, nopass :: calc_h1=>calc_hz
    procedure, nopass :: calc_h2=>calc_hx
    procedure, nopass :: calc_h3=>calc_hy
    procedure :: native2ECEFspher=>cart2ECEFspher

    !> type deallocations, reset flags, etc.
    final :: destructor
end type cartmesh


contains


!> allocate space and associate pointers with arrays in base class.  must runs set_coords first.
subroutine init_cartmesh(self)
  class(cartmesh), intent(inout) :: self

  if (.not. self%xi_alloc_status) error stop ' must have curvilinear coordinates defined prior to call init_cartmesh()'

  ! allocate array space using base type-bound procedure
  call self%calc_coord_diffs()
  print *,"cart:calc_coord_diffs done"
  call self%init_storage()
  print *, "cart:init_storage done"
  ! fixme: need to add geographic coord arrays first...
  !call self%calc_inull()

  ! now we must associate pointers for extended type alias variables.  This is mostly done in order
  !  to have more readable code.
  self%z=>self%x1; self%x=>self%x2; self%y=>self%x3
  self%zint=>self%x1i; self%xint=>self%x2i; self%yint=>self%x3i
  self%hz=>self%h1; self%hx=>self%h2; self%hy=>self%h3
  self%hzzi=>self%h1x1i; self%hxzi=>self%h2x1i; self%hyzi=>self%h3x1i
  self%hzxi=>self%h1x2i; self%hxxi=>self%h2x2i; self%hyxi=>self%h3x2i
  self%hzyi=>self%h1x3i; self%hxyi=>self%h2x3i; self%hyyi=>self%h3x3i
  self%ez=>self%e1; self%ex=>self%e2; self%ey=>self%e3
  self%gz=>self%g1; self%gx=>self%g2; self%gy=>self%g3
end subroutine init_cartmesh


!> create a cart mesh structure out of given q,p,phi spacings.  We assume here that the input cell center locations
!   are provide with ghost cells included (note input array indexing in dummy variable declarations.  For new we assume
!   that the fortran code will precompute and store the "full" grid information to save time (but this uses more memory).
subroutine make_cartmesh(self)
  class(cartmesh), intent(inout) :: self
  integer :: lzg,lxg,lyg,lz,lx,ly
  integer :: iz,ix,iy
  real(wp), dimension(:,:,:), pointer :: r,theta,phispher     ! so these can serve as targets, as needed by metric factor calculations
  real(wp) :: phictr,thetactr
  real(wp) :: gamma2,gamma1

  ! check that pointers are correctly associated, which implies that all space has been allocated :)
  if (.not. associated(self%z)) error stop  &
             ' pointers to grid coordiante arrays must be associated prior calling make_cartmesh()'

  ! size of arrays, including ghost cells
  lzg = size(self%z, 1)
  if(lzg < 1) error stop "lzg must be strictly positive"
  lxg = size(self%x, 1)
  if(lxg < 1) error stop "lxg must be strictly positive"
  lyg = size(self%y, 1)
  if(lyg < 1) error stop "lyg must be strictly positive"
  allocate(r(-1:lzg-2,-1:lxg-2,-1:lyg-2), theta(-1:lzg-2,-1:lxg-2,-1:lyg-2))
  allocate(phispher(-1:lzg-2,-1:lxg-2,-1:lyg-2))

  ! array sizes without ghost cells for convenience
  print '(A,1X,I0,1X,I0,1X,I0)', ' make_cartmesh:  allocating space for grid of size:  ',lzg,lxg,lyg
  lz=lzg-4; lx=lxg-4; ly=lyg-4;

  call self%native2ECEFspher(self%glonctr,self%glatctr,self%z,self%x,self%y,r,theta,phispher)

  ! now assign structure elements and deallocate unneeded temp variables
!  self%r=r(1:lz,1:lx,1:ly); self%theta=theta(1:lz,1:lx,1:ly); self%phi=phispher(1:lz,1:lx,1:ly)   ! don't need ghost cells!
  self%r=r(-1:lz+2,-1:lx+2,-1:ly+2); self%theta=theta(-1:lz+2,-1:lx+2,-1:ly+2); 
  self%phi=phispher(-1:lz+2,-1:lx+2,-1:ly+2)
  deallocate(r,theta,phispher)

  ! compute the geographic coordinates
  print*, ' make_cartmesh:  geographic coordinates from magnetic...'
  call self%calc_geographic()

  ! compute and store the metric factors; these need to include ghost cells
  !  these are function calls because I have to bind the deferred procedures...
  print*, ' make_cartmesh:  metric factors for cell centers...'
  self%hz(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h1(r,theta,phispher)
  self%hx(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h2(r,theta,phispher)
  self%hy(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h3(r,theta,phispher)

  ! q cell interface metric factors
  print*, ' make_cartmesh:  metric factors for cell q-interfaces...'
  self%hzzi=1
  self%hxzi=1
  self%hyzi=1

  ! p cell interface metric factors
  print*, ' make_cartmesh:  metric factors for cell p-intefaces...'
  self%hzxi=1
  self%hxxi=1
  self%hyxi=1

  print*, ' make_cartmesh:  metric factors for cell phi-interfaces...'
  self%hzyi=1
  self%hxyi=1
  self%hyyi=1

  ! spherical ECEF unit vectors (expressed in a Cartesian ECEF basis)
  print*, ' make_cartmesh:  spherical ECEF unit vectors...'
  call self%calc_er()
  call self%calc_etheta()
  call self%calc_ephi()

  ! cart coordinate system unit vectors (Cart. ECEF)
  print*, ' make_cartmesh:  cartesian unit vectors...'
  call self%calc_e1()
  call self%calc_e2()
  call self%calc_e3()

  ! magnetic field magnitude
  print*, ' make_cartmesh:  magnetic fields...'
  call self%calc_Bmag()

  ! gravity components
  print*, ' make_cartmesh:  gravity...'
  call self%calc_grav()

  ! set the status now that coord. specific calculations are done
  self%coord_alloc_status=.true.

  ! now finish by calling procedures from base abstract type
  print*, ' make_cartmesh:  base type-bound procedure calls...'
  call self%calc_difflengths()     ! differential lengths (units of m)
  call self%calc_inull()           ! null points (non computational)
  call self%calc_gridflag()        ! compute and store grid type

  ! inclination angle for each field line; awkwardly this must go after gridflag is set...
  print*, ' make_cartmesh:  inclination angle...'
  call self%calc_inclination()
  print *, "make_cartmesh done"
end subroutine make_cartmesh


!> utility convert native Cartesian coordinates to ECEF spherical
subroutine cart2ECEFspher(self,glonctr,glatctr,coord1,coord2,coord3,r,theta,phispher)
  class(cartmesh), intent(in) :: self
  real(wp) :: glonctr,glatctr
  real(wp), dimension(:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord2,1):ubound(coord2,1),lbound(coord3,1):ubound(coord3,1)), &
              intent(inout) :: r,theta,phispher
  real(wp), dimension(:), pointer :: x,y,z
  integer :: lx,ly,lz,ix,iy,iz
  integer :: izmin,izmax,ixmin,ixmax,iymin,iymax
  real(wp) :: thetactr,phictr
  real(wp) :: gamma1,gamma2

  z=>coord1; x=>coord2; y=>coord3;
  lz=size(z); lx=size(x); ly=size(y);

  izmin=lbound(z,1); izmax=ubound(z,1);   ! apparently pointers carry lbound/ubound into procedures but arrays don't???
  ixmin=lbound(x,1); ixmax=ubound(x,1);
  iymin=lbound(y,1); iymax=ubound(y,1);

  if (lz /= size(r,1) .or. lx /= size(r,2) .or. ly /= size(r,3) ) then
    print*, lz,lx,ly,size(r,1),size(r,2),size(r,3)
    error stop 'cartmesh::cart2ECEFspher - r and native coord. arrays not conformable...'
  end if

  call geog2geomag(glonctr,glatctr,phictr,thetactr)

  ! radial distance from Earth's center
  do iy=iymin,iymax
    do ix=ixmin,ixmax
      do iz=izmin,izmax
        r(iz,ix,iy)=Re+z(iz)
      end do
    end do
  end do

  ! northward angular distance
  do iy=iymin,iymax
    gamma2=y(iy)/Re                  ! must retain sign(y)
    do ix=ixmin,ixmax
      do iz=izmin,izmax
        theta(iz,ix,iy)=thetactr-gamma2  ! minus because theta is positive south (y runs positive north)
      end do
    end do
  end do

  ! eastward angular distance
  do iy=iymin,iymax
    do ix=ixmin,ixmax
      gamma1=x(ix)/Re/sin(thetactr)
      do iz=izmin,izmax
        phispher(iz,ix,iy)=phictr+gamma1
      end do
    end do
  end do
end subroutine cart2ECEFspher


!> compute gravitational field components
subroutine calc_grav_cart(self)
  class(cartmesh), intent(inout) :: self
  ! fixme: error checking?

  print*, size(self%r,1),size(self%r,2),size(self%r,3)
!  self%gz=-Gconst*Me/self%r(1:size(self%r,1)-4,1:size(self%r,2)-4,1:size(self%r,3)-4)**2     ! radial component of gravity
  self%gz=-Gconst*Me/self%r(1:self%lx1,1:self%lx2,1:self%lx3)**2     ! radial component of gravity
  self%gx=0
  self%gy=0
end subroutine calc_grav_cart


!> compute the magnetic field strength.  For cartesian grids assume constant
subroutine calc_Bmag_cart(self)
  class(cartmesh), intent(inout) :: self

  ! fixme: error checking

  !self%Bmag=mu0*Mmag/4/pi/self%r**3*sqrt(3*cos(self%theta)**2+1)
  self%Bmag=-50000.e-9_wp
end subroutine calc_Bmag_cart


!> compute the inclination angle (degrees) for each geomagnetic field line
!   for a Cartesian grid the only thing that really makes any sense is 90 degrees
subroutine calc_inclination_cart(self)
  class(cartmesh), intent(inout) :: self

  ! fixme: error checking

  self%I=90
end subroutine calc_inclination_cart


!> compute metric factors for q
function calc_hz(coord1,coord2,coord3) result(hval)
  !! FIXME: add error checking
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi


  logical :: a
  a = associated(coord2)
  !! avoid unused variable warning

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1
end function calc_hz


!> compute p metric factors
function calc_hx(coord1,coord2,coord3) result(hval)
  !! FIXME: add error checking
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  logical :: a
  a = associated(coord2)
  !! avoid unused variable warning

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1
end function calc_hx


!> compute phi metric factor
function calc_hy(coord1,coord2,coord3) result(hval)
  !! FIXME: add error checking
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  logical :: a
  a = associated(coord2)
  !! avoid unused variable warning

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1
end function calc_hy


!> radial unit vector (expressed in ECEF cartesian coodinates, components permuted as ix,iy,iz)
subroutine calc_er_spher(self)
  class(cartmesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%er=er_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_er_spher


!> zenith angle unit vector (expressed in ECEF cartesian coodinates
subroutine calc_etheta_spher(self)
  class(cartmesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%etheta=etheta_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_etheta_spher


!> azimuth angle unit vector (ECEF cart.)
subroutine calc_ephi_spher(self)
  class(cartmesh), intent(inout) :: self
  integer :: lx1,lx2,lx3

  lx1=self%lx1; lx2=self%lx2; lx3=self%lx3;

  ! fixme: error checking

  self%ephi=ephi_spherical(self%theta(1:lx1,1:lx2,1:lx3),self%phi(1:lx1,1:lx2,1:lx3))
end subroutine calc_ephi_spher


!> unit vector in the z direction
subroutine calc_ez(self)
  class(cartmesh), intent(inout) :: self

  ! fixme: error checking

  self%ez=self%er     ! altitude is same direction as distance from Earth's center
end subroutine calc_ez


!> unit vector in the p direction
subroutine calc_ex(self)
  class(cartmesh), intent(inout) :: self

  ! fixme: error checking

  self%ex=self%ephi    ! x-direction is eastward
end subroutine calc_ex


!> unit vector in the phi direction
subroutine calc_ey(self)
  class(cartmesh), intent(inout) :: self

  ! fixme: error checking

  self%ey=-self%etheta    ! ey should be positive north (etheta is positive south)
end subroutine calc_ey


!> type destructor; written generally, viz. as if it is possible some grid pieces are allocated an others are not
subroutine destructor(self)
  type(cartmesh) :: self

  call self%dissociate_pointers()
  print*, '  cartmesh destructor completed successfully'
end subroutine destructor

end module meshobj_cart
