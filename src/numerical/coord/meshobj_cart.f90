module meshobj_cart

!> Contains data and subroutines for managing a Cartesian mesh

! uses
use phys_consts, only: wp,Re,pi,Mmag,mu0,Gconst,Me
use meshobj, only: curvmesh
use spherical, only: er_spherical,etheta_spherical,ephi_spherical
use geomagnetic, only: geog2geomag,geomag2geog,r2alt,alt2r

implicit none (type, external)


! type extension for dipolemesh
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

  real(wp) :: glonctr,glatctr     ! center of the grid in geographic lat/lon

  contains
    !> type-bound procs. for dipole meshes
    procedure :: set_center
 
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
    procedure :: calc_geographic=>calc_geographic_cart
    
    !> type deallocations, reset flags, etc.
    final :: destructor
end type cartmesh


contains


!> allocate space and associate pointers with arrays in base class.  must runs set_coords first.
subroutine init_cartmesh(self)
  class(cartmesh), intent(inout) :: self

  if (.not. self%xi_alloc_status) error stop ' must have curvilinear coordinates defined prior to call init_dipolemesh()'

  ! allocate array space using base type-bound procedure
  call self%calc_coord_diffs()
  call self%init_storage()
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


!> create a dipole mesh structure out of given q,p,phi spacings.  We assume here that the input cell center locations
!   are provide with ghost cells included (note input array indexing in dummy variable declarations.  For new we assume
!   that the fortran code will precompute and store the "full" grid information to save time (but this uses more memory).  
subroutine make_cartmesh(self) 
  class(cartmesh), intent(inout) :: self
  integer :: lzg,lxg,lyg,lz,lx,ly
  integer :: iz,ix,iy
  real(wp), dimension(:,:,:), pointer :: r,theta,phispher     ! so these can serve as targets, as needed by metric factor calculations
  real(wp) :: gamma2,gamma1

  ! check that pointers are correctly associated, which implies that all space has been allocated :)
  if (.not. associated(self%z)) error stop  & 
             ' pointers to grid coordiante arrays must be associated prior calling make_dipolemesh()'

  ! size of arrays, including ghost cells
  lzg=size(self%z,1); lxg=size(self%x,1); lyg=size(self%y,1)
  allocate(r(-1:lzg-2,-1:lxg-2,-1:lyg-2),theta(-1:lzg-2,-1:lxg-2,-1:lyg-2))
  allocate(phispher(-1:lzg-2,-1:lxg-2,-1:lyg-2))
  
  ! array sizes without ghost cells for convenience
  print*, ' make_dipolemesh:  allocating space for grid of size:  ',lzg,lxg,lyg
  lz=lzg-4; lx=lxg-4; ly=lyg-4;

  ! convert the cell centers to spherical ECEF coordinates, then tile for longitude dimension
  print*, ' make_dipolemesh:  converting cell centers to spherical coordinates...'
  call geog2geomag(self%glonctr,self%glatctr,phictr,thetactr)

  ! radial distance from Earth's center
  do iy=-1,ly+2
    do ix=-1,lx+2
      do iz=-1,lz+2
        r(iz,ix,iy)=Re+z(iz)
      end do
    end do
  end do

  ! northward angular distance
  do iy=-1,ly+2
    gamma2=y(iy)/Re                  ! must retain sign(y)
    do ix=-1,lx+2
      do iz=-1,lz+2
        theta(iz,ix,iy)=thetactr-gamma2  ! minus because theta is positive south (y runs positive north)
      end do
    end do
  end do

  ! eastward angular distance
  do iy=-1,ly+2
    do ix=-1,lx+2
      gamma1=x(ix)/Re/sin(thetactr)
      do iz=-1,lz+2
        phispher(iz,ix,iy)=phictr+gamma1
      end do
    end do
  end do

  ! now assign structure elements and deallocate unneeded temp variables
  self%r=r(1:lz,1:lx,1:ly); self%theta=theta(1:lz,1:lx,1:ly); self%phi=phispher(1:lx,1:lx,1:ly)   ! don't need ghost cells!
  deallocate(r,theta,phispher)

  ! compute the geographic coordinates
  print*, ' make_dipolemesh:  geographic coordinates from magnetic...'
  call self%calc_geographic() 

  ! compute and store the metric factors; these need to include ghost cells
  !  these are function calls because I have to bind the deferred procedures...
  print*, ' make_dipolemesh:  metric factors for cell centers...'
  self%hz(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h1(r,theta,phispher)
  self%hx(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h2(r,theta,phispher)
  self%hy(-1:lz+2,-1:lx+2,-1:ly+2)=self%calc_h3(r,theta,phispher)

  ! q cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell q-interfaces...'
  self%hzzi=1._wp
  self%hxzi=1._wp
  self%hyzi=1._wp

  ! p cell interface metric factors
  print*, ' make_dipolemesh:  metric factors for cell p-intefaces...'
  self%hzxi=1._wp
  self%hxxi=1._wp
  self%hyxi=1._wp

  print*, ' make_dipolemesh:  metric factors for cell phi-interfaces...'
  self%hzyi=1._wp
  self%hxyi=1._wp
  self%hyyi=1._wp

  ! spherical ECEF unit vectors (expressed in a Cartesian ECEF basis)
  print*, ' make_dipolemesh:  spherical ECEF unit vectors...'  
  call self%calc_er()
  call self%calc_etheta()
  call self%calc_ephi()

  ! dipole coordinate system unit vectors (Cart. ECEF)
  print*, ' make_dipolemesh:  dipole unit vectors...'  
  call self%calc_e1()
  call self%calc_e2()
  call self%calc_e3()

  ! magnetic field magnitude
  print*, ' make_dipolemesh:  magnetic fields...'    
  call self%calc_Bmag()

  ! gravity components
  print*, ' make_dipolemesh:  gravity...'
  call self%calc_grav()

  ! set the status now that coord. specific calculations are done
  self%coord_alloc_status=.true.

  ! now finish by calling procedures from base abstract type
  print*, ' make_dipolemesh:  base type-bound procedure calls...'
  call self%calc_difflengths()     ! differential lengths (units of m)
  call self%calc_inull()           ! null points (non computational)
  call self%calc_gridflag()        ! compute and store grid type

  ! inclination angle for each field line; awkwardly this must go after gridflag is set...
  print*, ' make_dipolemesh:  inclination angle...'  
  call self%calc_inclination()
end subroutine make_cartmesh


!> this sets the geographic center of the mesh; must be called prior to make()
subroutine set_center(self,glon,glat)
  class(cartmesh), intent(inout) :: self
  real(wp), intent(in) :: glon,glat

  self%glonctr=glon
  self%glatctr=glat
end subroutine set_center


!> compute geographic coordinates of all grid points
subroutine calc_geographic_cart(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking?

  call geomag2geog(self%phi,self%theta,self%glon,self%glat)
  self%alt=r2alt(self%r)
  self%geog_set_status=.true.
end subroutine calc_geographic_cart


!> compute gravitational field components
subroutine calc_grav_cart(self)
  class(dipolemesh), intent(inout) :: self
  real(wp), dimension(1:self%lx1,1:self%lx2,1:self%lx3) :: gr
 
  ! fixme: error checking?
 
  self%gz=-Gconst*Me/self%r**2     ! radial component of gravity
  self%gx=0._wp
  self%gy=0._wp
end subroutine calc_grav_cart


!> compute the magnetic field strength
subroutine calc_Bmag_cart(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%Bmag=mu0*Mmag/4/pi/self%r**3*sqrt(3*cos(self%theta)**2+1)
end subroutine calc_Bmag_cart


!> compute the inclination angle (degrees) for each geomagnetic field line
!   for a Cartesian grid the only thing that really makes any sense is 90 degrees
subroutine calc_inclination_dipole(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%I=90._wp
end subroutine calc_inclination_dipole


!> compute metric factors for q 
function calc_hz(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checking

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1._wp
end function calc_hz


!> compute p metric factors
function calc_hx(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checkign

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1._wp
end function calc_hx


!> compute phi metric factor
function calc_hy(coord1,coord2,coord3) result(hval)
  real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
  real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                      lbound(coord1,3):ubound(coord1,3)) :: hval
  real(wp), dimension(:,:,:), pointer :: r,theta,phi

  ! fixme: error checking

  r=>coord1; theta=>coord1; phi=>coord3;
  hval=1._wp
end function calc_hy


!> radial unit vector (expressed in ECEF cartesian coodinates, components permuted as ix,iy,iz)
subroutine calc_er_spher(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%er=er_spherical(self%theta,self%phi)
end subroutine calc_er_spher


!> zenith angle unit vector (expressed in ECEF cartesian coodinates
subroutine calc_etheta_spher(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%etheta=etheta_spherical(self%theta,self%phi)
end subroutine calc_etheta_spher


!> azimuth angle unit vector (ECEF cart.)
subroutine calc_ephi_spher(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%ephi=ephi_spherical(self%theta,self%phi)
end subroutine calc_ephi_spher


!> unit vector in the z direction
subroutine calc_ez(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%ez=self%er     ! altitude is same direction as distance from Earth's center
end subroutine calc_ez


!> unit vector in the p direction
subroutine calc_ex(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%ex=self%ephi    ! x-direction is eastward
end subroutine calc_ex


!> unit vector in the phi direction
subroutine calc_ey(self)
  class(dipolemesh), intent(inout) :: self

  ! fixme: error checking

  self%ey=-self%etheta    ! ey should be positive north (etheta is positive south)
end subroutine calc_ey


!> type destructor; written generally, viz. as if it is possible some grid pieces are allocated an others are not
subroutine destructor(self)
  type(dipolemesh) :: self

  call self%dissociate_pointers()
  print*, '  cartmesh destructor completed successfully'
end subroutine destructor

end module meshobj_cart

