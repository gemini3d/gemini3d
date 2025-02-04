module meshobj

!> This is a base type for defining most general characteristics of a curvlinear mesh; extension define specific coordinate systems.
!    However, the idea here is to insulate the numerical parts of the program from that so that they can deal exclusively with the
!    generic metric factors, etc.

use phys_consts, only : wp,pi
use h5fortran, only : hdf5_file
use geomagnetic, only: geog2geomag,geomag2geog,r2alt,alt2r,rotgg2gm
use spherical, only: er_spherical,etheta_spherical,ephi_spherical

implicit none (type, external)
private
public :: curvmesh

!> curvmesh is an abstract type containing functionality and data that is not specific to individual coordinate systems
!   (which are extended types).  Note that all arrays are pointers because they need to be targets and allocatable AND the fortran
!   standard does not support having allocatable, target attributed inside a derived type.  Because of this is it not straightforward
!   to check the allocation status of these arrays (i.e. fortran also does not allow one to check
!   the allocation status of a pointer).  Thus the quantities which are not, for sure, allocated need to have an allocation status
!   variable so we can check...  Because the pointers are always allocated in groups we do not need separate status vars for each array
!   thankfully...
!
!  Note that the deferred bindings all have the single argument of self so any data used to compute whatever quantity is desired must
!   be stored in the derived type.  This is done to avoid long, potentially superfluous argument lists that may contain extraneous information,
!   e.g. if a coordinate is not needed to compute a metric factor and similar situations.  Any
!   operation that needs more input should be defined as a type-bound procedure in any extensions to this abstract type.  This also means that
!   error checking needs to be done in the extended type to insure that the data needed for a given calculation exist.
type, abstract :: curvmesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Generic properties !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> we need to know whether or not different groups of pointers are allocated and the intrinsic
  !  "allocatable" only work on allocatables (not pointers)...
  logical :: xi_alloc_status=.false.
  logical :: dxi_alloc_status=.false.
  logical :: dxi_alloc_status_root=.false.
  logical :: difflen_alloc_status=.false.
  logical :: null_alloc_status=.false.
  logical :: geog_set_status=.false.      ! geographic coords. get allocated with other arrays, but get set separately
  logical :: coord_set_status_root=.false.
  logical :: geogi_set_status=.false.

  !> sizes.  Specified and set by base class methods
  integer :: lx1,lx2,lx3,lx2all,lx3all
  !! for program units that may not be able to access module globals

  !!> curvilinear vars. and diffs.
  real(wp), dimension(:), pointer :: x1     ! provided in input file
  real(wp), dimension(:), pointer :: x1i    ! recomputed by base class once x1 set
  real(wp), dimension(:), pointer :: dx1        ! recomputed by base class from x1
  !! because sub arrays need to be assigned to aliases in calculus module program units
  real(wp), dimension(:), pointer :: dx1i   ! recomputed by base class from interface locations (themselves recomputed)

  real(wp), dimension(:), pointer :: x2
  real(wp), dimension(:), pointer :: x2i
  real(wp), dimension(:), pointer :: dx2
  !! this (and similar dx arrays) are pointers because they become associated with
  !! other pointers (must either have the "target" keyword or are themselves pointers).
  !! These should also be contiguous but I believe that is guaranteed as long as they are assigned
  !! through allocate statements (see fortran standard)  derived type arrays cannot be targets so we
  !! are forced to declare them as pointers and trust that they are allocated contiguous
  real(wp), dimension(:), pointer :: dx2i

  real(wp), dimension(:), pointer :: x3
  real(wp), dimension(:), pointer :: x3i
  real(wp), dimension(:), pointer :: dx3
  real(wp), dimension(:), pointer :: dx3i

  !> these are fullgrid but possibly carried around by workers too since 1D arrays?
  real(wp), dimension(:), pointer :: x2all
  real(wp), dimension(:), pointer  :: x2iall
  real(wp), dimension(:), pointer :: dx2all
  real(wp), dimension(:), pointer  :: dx2iall

  real(wp), dimension(:), pointer :: x3all
  real(wp), dimension(:), pointer  :: x3iall
  real(wp), dimension(:), pointer :: dx3all
  real(wp), dimension(:), pointer  :: dx3iall

  !> differential length elements.  Compute by a method in generic class given metric factors
  real(wp), dimension(:,:,:), pointer :: dl1i,dl2i,dl3i

  !> flag for indicating whether or not the grid is periodic
  logical :: flagper=.false.

  !> flag for indicated type of grid (0 - closed dipole; 1 - open dipole inverted; 2 - non-inverted).
  !! Computed by method in abstract type once coordinate specific quantitaties are computed by extension.
  integer :: gridflag

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Coordinate system specific properties !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Metric factors.  These are pointers to be assigned/allocated/filled by type extensions for
  !! specific coordinate system
  logical :: coord_alloc_status=.false.    ! single status variable for all coord-specific arrays
  logical :: coord_alloc_status_root=.false.

  real(wp), dimension(:,:,:), pointer :: h1,h2,h3     ! need to be computed by subclass for specific coordinate system
  !! these are the cell-centered metric coefficients
  real(wp), dimension(:,:,:), pointer :: h1x1i,h2x1i,h3x1i
  !! metric factors at x1 cell interfaces; dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), pointer :: h1x2i,h2x2i,h3x2i
  !! metric factors at x2 interfaces; dim. 2 has lx2+1
  real(wp), dimension(:,:,:), pointer :: h1x3i,h2x3i,h3x3i
  !! metric factors at x3 interfaces; dim. 3 has lx3+1

  !> root only; full grid metric factors.  These must have type-extension bound procedures to
  !> provide allocation and filling
  real(wp), dimension(:,:,:), pointer :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), pointer :: h1x1iall,h2x1iall,h3x1iall
  !! dimension 1 has size lx1+1
  real(wp), dimension(:,:,:), pointer :: h1x2iall,h2x2iall,h3x2iall
  !! dim. 2 has lx2all+1
  real(wp), dimension(:,:,:), pointer :: h1x3iall,h2x3iall,h3x3iall
  !! dim. 3 has lx3all+1

  !> Problem-depedent geometric terms that can be precomputed, e.g. for advection and elliptic
  !! equations
  !! to be added in later if we deem it worthwhile to include these to save computation time and/or
  !! memory

  !> unit vectors - Pointers.  Subclass methods must allocate and assign values to these.
  real(wp), dimension(:,:,:,:), pointer :: e1,e2,e3
  !! unit vectors in curvilinear space (cartesian components)
  real(wp), dimension(:,:,:,:), pointer :: er,etheta,ephi
  !! spherical unit vectors (cartesian components)

  !> geomagnetic grid data, ECEF spherical referenced to dipole axis
  real(wp), dimension(:,:,:), pointer :: r,theta,phi
  !! may be used in the interpolation of neutral perturbations

  !> root only; full grid geomagnetic positions.  Pointers.  Need type-extension bound procedures
  !! for each
  !coordinate system must allocate and compute these.
  real(wp), dimension(:,:,:), pointer :: rall,thetaall,phiall

  !> geographic data; pointers
  real(wp), dimension(:,:,:), pointer :: glat,glon,alt
  real(wp), dimension(:,:,:), pointer :: glati,gloni,alti

  !> magnetic field magnitude and inclination.  Pointers
  real(wp), dimension(:,:,:), pointer :: Bmag
  real(wp), dimension(:,:), pointer :: I

  !> Gravitational field
  real(wp), dimension(:,:,:), pointer :: g1,g2,g3

  !> EIA-required fullgrid data.  Pointers
  real(wp), dimension(:,:,:), pointer :: altall,glonall
  real(wp), dimension(:,:,:), pointer :: Bmagall

  !> points that should not be considered part of the numerical domain. Pointers
  logical, dimension(:,:,:), pointer :: nullpts
  !! this could be a logical but I'm going to treat it as real*8
  integer :: lnull
  !! length of null point index array
  integer, dimension(:,:), pointer :: inull

  !> for floating coordinate systems we must define an aboslute center position in geographic
  !! coordinates
  real(wp) :: glonctr,glatctr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! type-bound procedures !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  contains
    procedure :: set_coords             ! initialize general curvilinear coordinates and mesh sizes
    procedure :: calc_coord_diffs       ! compute bwd and midpoint diffs from coordinates
    procedure :: calc_coord_diffs_root  ! coordinate diffs for root fullgrid
    procedure :: calc_difflengths       ! compute differential lengths
    procedure :: calc_inull             ! compute null points
    procedure :: calc_gridflag          ! compute the type of grid we have
    procedure :: init_storage           ! allocate space for coordinate specific arrays
    procedure :: init_storage_root
    procedure :: writesize              ! write a simsize.h5 file for this local grid
    procedure :: writegrid              ! write curvilinear coordinates to simgrid.h5
    procedure :: writegridall           ! write all mesh arrays to simgrid.h5
    procedure :: dissociate_pointers    ! clear out memory and reset allocation flags
    procedure :: calc_geographic        ! convert to geographic coordinates
    procedure :: set_root               ! set fullgrid variables that have been gathered from workers to root
    procedure :: set_periodic           ! set the flag which labels grid as periodic vs. aperiodic
    procedure :: set_center             ! set the center of the grid (only needed for "floating" coordinate systems)
    procedure :: calc_unitvec_geo       ! compute geographic unit vectors over the grid
    procedure :: calc_geographici       ! cell edge locations
    !final :: destructor
    !! an abstract type cannot have a final procedure, as the final procedure
    !! must act on a type and not polymorhpic object

    !! deferred bindings and associated generic interfaces
    procedure(initmake), deferred :: init
    procedure(initmake), deferred :: make
    !! There is no make_root procedure as this requires message passing and needs to be done outside
    !! this type
    procedure(calc_procedure), deferred :: calc_grav
    procedure(calc_procedure), deferred :: calc_Bmag
    procedure(calc_procedure), deferred :: calc_inclination
    procedure(calc_procedure), deferred :: calc_er
    procedure(calc_procedure), deferred :: calc_etheta
    procedure(calc_procedure), deferred :: calc_ephi
    procedure(calc_procedure), deferred :: calc_e1
    procedure(calc_procedure), deferred :: calc_e2
    procedure(calc_procedure), deferred :: calc_e3
    procedure(native_convert), deferred :: native2ECEFspher

    !! these bindings have different interfaces due to evaluation of metric factors at cell centers
    !! vs. interfaces
    !   because the return value may need to be assigned to different member arrays this should be a
    !!   function that
    !   does not pass in or operate on its object.
    procedure(calc_metric), deferred, nopass :: calc_h1
    procedure(calc_metric), deferred, nopass :: calc_h2
    procedure(calc_metric), deferred, nopass :: calc_h3
end type curvmesh


!> interfaces for deferred bindings, most will operate on data in self - this provides maximum
!! flexibility for the extension to use whatever data it needs to compute the various grid quantities
abstract interface
  subroutine initmake(self)
    import curvmesh
    class(curvmesh), intent(inout) :: self
  end subroutine initmake
  subroutine calc_procedure(self)
    import curvmesh
    class(curvmesh), intent(inout) :: self
  end subroutine calc_procedure
  function calc_metric(coord1,coord2,coord3) result(hval)
    import wp
    real(wp), dimension(:,:,:), pointer, intent(in) :: coord1,coord2,coord3
    !!! because these will be aliased (readability) so should be pointers, e.g. x,y,z or r,theta,phi
    real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
                        lbound(coord1,3):ubound(coord1,3)) :: hval   ! arrays may not start at index 1
  end function calc_metric
  subroutine native_convert(self,glonctr,glatctr,coord1,coord2,coord3,r,theta,phispher)
    import curvmesh
    import wp
    class(curvmesh), intent(in) :: self
    real(wp) :: glonctr,glatctr
    real(wp), dimension(:), pointer, intent(in) :: coord1,coord2,coord3
    !real(wp), dimension(:,:,:), intent(inout) :: r,theta,phispher
    real(wp), dimension(lbound(coord1,1):ubound(coord1,1),lbound(coord2,1):ubound(coord2,1), &
                lbound(coord3,1):ubound(coord3,1)), &
                intent(inout) :: r,theta,phispher
  end subroutine native_convert
end interface


contains
  subroutine set_coords(self,x1,x2,x3,x2all,x3all)
  !! assign coordinates to internal variables given some set of input arrays.
  !!  Assume that the data passed in include ghost cells
    class(curvmesh), intent(inout) :: self
    real(wp), dimension(:), intent(in) :: x1,x2,x3
    real(wp), dimension(:), intent(in) :: x2all,x3all

    integer :: lx1,lx2,lx3,lx2all,lx3all

    lx1=size(x1,1)-4
    lx2=size(x2,1)-4
    lx3=size(x3,1)-4
    lx2all=size(x2all,1)-4
    lx3all=size(x3all,1)-4

    if(lx1 < 1) error stop 'meshobj:set_coords: lx1 must be strictly positive'
    if(lx2 < 1) error stop 'meshobj:set_coords: lx2 must be strictly positive'
    if(lx3 < 1) error stop 'meshobj:set_coords: lx3 must be strictly positive'
    if(lx2all < lx2) error stop 'meshobj:set_coords: lx2all must be > lx2'
    if(lx3all < lx3) error stop 'meshobj:set_coords: lx3all must be > lx3'

    self%lx1=lx1
    self%lx2=lx2
    self%lx3=lx3
    self%lx2all=lx2all
    self%lx3all=lx3all

    allocate(self%x1(-1:lx1+2), self%x2(-1:lx2+2), self%x3(-1:lx3+2))
    self%x1=x1
    self%x2=x2
    self%x3=x3
    allocate(self%x2all(-1:lx2all+2), self%x3all(-1:lx3all+2))
    self%x2all=x2all
    self%x3all=x3all

    self%xi_alloc_status=.true.
  end subroutine set_coords


  subroutine set_root(self,h1all,h2all,h3all, &
                      h1x1iall,h2x1iall,h3x1iall, &
                      h1x2iall,h2x2iall,h3x2iall, &
                      h1x3iall,h2x3iall,h3x3iall, &
                      rall,thetaall,phiall, &
                      altall,Bmagall,glonall)
  !! Set full grid arrays needed for some numerical parts - these need to be gathered by
  !!   a wrapping module/procedure and then passed to this procedure.
    class(curvmesh), intent(inout) :: self
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: h1all,h2all,h3all
    real(wp), dimension(:,:,:), intent(in) :: h1x1iall,h2x1iall,h3x1iall
    real(wp), dimension(:,:,:), intent(in) :: h1x2iall,h2x2iall,h3x2iall
    real(wp), dimension(:,:,:), intent(in) :: h1x3iall,h2x3iall,h3x3iall
    real(wp), dimension(:,:,:), intent(in) :: rall,thetaall,phiall
    real(wp), dimension(:,:,:), intent(in) :: altall,Bmagall,glonall

    if (.not. self%coord_alloc_status_root) error stop ' attempting to set root params. before allocating!'

    self%h1all=h1all
    self%h2all=h2all
    self%h3all=h3all
    self%h1x1iall=h1x1iall
    self%h2x1iall=h2x1iall
    self%h3x1iall=h3x1iall
    self%h1x2iall=h1x2iall
    self%h2x2iall=h2x2iall
    self%h3x2iall=h3x2iall
    self%h1x3iall=h1x3iall
    self%h2x3iall=h2x3iall
    self%h3x3iall=h3x3iall
    self%rall=rall
    self%thetaall=thetaall
    self%phiall=phiall
    self%altall=altall
    self%Bmagall=Bmagall
    self%glonall=glonall
    self%coord_set_status_root=.true.
  end subroutine set_root


  !> assign a reference meridian alt,lat,lon across periodic grid
  subroutine set_periodic(self,flagperiodic,refalt,refglon,refglat)
    class(curvmesh), intent(inout) :: self
    integer, intent(in) :: flagperiodic
    real(wp), dimension(:,:), intent(in) :: refalt,refglon,refglat
    integer :: ix3

    ! flag appropriately
    if (flagperiodic/=0) then
      self%flagper=.true.
    else
      self%flagper=.false.
    end if

    ! In the special case where the user wants flagperiodic==1, we take the additional step of forcing
    !   glat/glon to be the same across the x3 dimension so that the neutral atmosphere and photoionization
    !   will be constant vs. x3.  This is most typically used in instability simulations when one wants to
    !   explicitly remove any dependence of background parameters on the x3 direction.  One would not want to use
    !   this when simply trying to model an angular coordinate (e.g. longitude) across the full globe as the
    !   atmospheric and SZA changes are needed to realistically capture the system.
    ! force periodicity in geographic locations using reference meridian data
    !! FIXME: specify ghost cells for glat/glon arrays?
    if (flagperiodic==1) then
      do ix3=1,self%lx3
        self%glat(1:self%lx1,1:self%lx2,ix3)=refglat(:,:)
        self%glon(1:self%lx1,1:self%lx2,ix3)=refglon(:,:)
        self%alt(1:self%lx1,1:self%lx2,ix3)=refalt(:,:)
      end do
    end if
  end subroutine set_periodic


  !> compute diffs from given grid spacing.  Note that no one except for root needs full grid diffs.
  subroutine calc_coord_diffs(self)
    class(curvmesh), intent(inout) :: self

    integer :: lx1,lx2,lx3

    if (.not. self%xi_alloc_status) error stop ' attempting to compute diffs without coordinates!'

    lx1=self%lx1
    lx2=self%lx2
    lx3=self%lx3

    if(lx1 < 1) error stop 'meshobj:calc_coord_diffs: lx1 must be strictly positive'
    if(lx2 < 1) error stop 'meshobj:calc_coord_diffs: lx2 must be strictly positive'
    if(lx3 < 1) error stop 'meshobj:calc_coord_diffs: lx3 must be strictly positive'

    allocate(self%dx1(0:lx1+2), self%x1i(1:lx1+1), self%dx1i(1:lx1))
    self%dx1 = self%x1(0:lx1+2)-self%x1(-1:lx1+1)
    !! computing these avoids extra message passing (could be done for other coordinates
    self%x1i(1:lx1+1) = 0.5_wp*(self%x1(0:lx1)+self%x1(1:lx1+1))
    self%dx1i=self%x1i(2:lx1+1)-self%x1i(1:lx1)

    allocate(self%dx2(0:lx2+2), self%x2i(1:lx2+1), self%dx2i(1:lx2))
    self%dx2 = self%x2(0:lx2+2)-self%x2(-1:lx2+1)
    !! computing these avoids extra message passing (could be done for other coordinates
    self%x2i(1:lx2+1) = 0.5_wp*(self%x2(0:lx2)+self%x2(1:lx2+1))
    self%dx2i=self%x2i(2:lx2+1)-self%x2i(1:lx2)

    allocate(self%dx3(0:lx3+2), self%x3i(1:lx3+1), self%dx3i(1:lx3))
    self%dx3 = self%x3(0:lx3+2)-self%x3(-1:lx3+1)
    self%x3i(1:lx3+1)=0.5_wp*(self%x3(0:lx3)+self%x3(1:lx3+1))
    self%dx3i=self%x3i(2:lx3+1)-self%x3i(1:lx3)

    self%dxi_alloc_status=.true.
  end subroutine calc_coord_diffs


  subroutine calc_coord_diffs_root(self)
    class(curvmesh), intent(inout) :: self

    integer :: lx1,lx2all,lx3all

    if (.not. self%xi_alloc_status) error stop ' attempting to compute root diffs without coordinates!'

    lx1=self%lx1
    lx2all=self%lx2all
    lx3all=self%lx3all

    if(lx1 < 1) error stop 'meshobj:calc_coord_diffs_root: lx1 must be strictly positive'
    if(lx2all < 1) error stop 'meshobj:calc_coord_diffs_root: lx2all must be strictly positive'
    if(lx3all < 1) error stop 'meshobj:calc_coord_diffs_root: lx3all must be strictly positive'

    allocate(self%x2iall(1:lx2all+1),self%dx2all(0:lx2all+2),self%dx2iall(1:lx2all))
    self%dx2all = self%x2all(0:lx2all+2)-self%x2all(-1:lx2all+1)
    self%x2iall(1:lx2all+1) = 0.5_wp*(self%x2all(0:lx2all)+self%x2all(1:lx2all+1))
    self%dx2iall=self%x2iall(2:lx2all+1)-self%x2iall(1:lx2all)
    allocate(self%x3iall(1:lx3all+1),self%dx3all(0:lx3all+2),self%dx3iall(1:lx3all))
    self%dx3all = self%x3all(0:lx3all+2)-self%x3all(-1:lx3all+1)
    self%x3iall(1:lx3all+1)=0.5_wp*(self%x3all(0:lx3all)+self%x3all(1:lx3all+1))
    self%dx3iall=self%x3iall(2:lx3all+1)-self%x3iall(1:lx3all)

    self%dxi_alloc_status_root=.true.
  end subroutine calc_coord_diffs_root


  !> calculate differential lengths (units of m), needed for CFL calculations
  subroutine calc_difflengths(self)
    class(curvmesh) :: self
    real(wp), dimension(:,:,:), allocatable :: tmpdx
    integer :: lx1,lx2,lx3

    if (.not. self%dxi_alloc_status .or. .not. self%coord_alloc_status) then
      print*, self%dxi_alloc_status,self%coord_alloc_status
      print*, (.not. self%dxi_alloc_status), (.not. self%coord_alloc_status)
      print*, .not. self%dxi_alloc_status .eqv. self%dxi_alloc_status
      error stop ' attempting to compute differential lengths without interface diffs or metric factors!'
    end if

    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    allocate(tmpdx(1:lx1,1:lx2,1:lx3))
    allocate(self%dl1i(1:lx1,1:lx2,1:lx3),self%dl2i(1:lx1,1:lx2,1:lx3),self%dl3i(1:lx1,1:lx2,1:lx3))
    tmpdx=spread(spread(self%dx1i,2,lx2),3,lx3)
    self%dl1i=tmpdx*self%h1(1:lx1,1:lx2,1:lx3)
    tmpdx=spread(spread(self%dx2i,1,lx1),3,lx3)
    self%dl2i=tmpdx*self%h2(1:lx1,1:lx2,1:lx3)
    tmpdx=spread(spread(self%dx3i,1,lx1),2,lx2)
    self%dl3i=tmpdx*self%h3(1:lx1,1:lx2,1:lx3)
    deallocate(tmpdx)

    self%difflen_alloc_status=.true.
  end subroutine calc_difflengths


  !> allocate space for metric factors, unit vectors, and transformations
  subroutine init_storage(self)
    class(curvmesh), intent(inout) :: self

    integer :: lx1,lx2,lx3

    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    if(lx1 < 1) error stop 'meshobj:init_storage: lx1 must be strictly positive'
    if(lx2 < 1) error stop 'meshobj:init_storage: lx2 must be strictly positive'
    if(lx3 < 1) error stop 'meshobj:init_storage: lx3 must be strictly positive'

    if (self%coord_alloc_status) error stop 'attempting to allocate space for coordinate-specific arrays when they already exist!'
    !! use this as a proxy for if any other coordinate-specific arrays exist
    !allocate(self%h1(1:lx1,1:lx2,1:lx3),self%h2(1:lx1,1:lx2,1:lx3),self%h3(1:lx1,1:lx2,1:lx3))
    allocate(self%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(self%h2, self%h3, mold=self%h1)

    allocate(self%h1x1i(1:lx1+1,1:lx2,1:lx3))
    allocate(self%h2x1i, self%h3x1i, mold=self%h1x1i)

    allocate(self%h1x2i(1:lx1,1:lx2+1,1:lx3))
    allocate(self%h2x2i, self%h3x2i, mold=self%h1x2i)

    allocate(self%h1x3i(1:lx1,1:lx2,1:lx3+1))
    allocate(self%h2x3i,self%h3x3i, mold=self%h1x3i)

    allocate(self%er(1:lx1,1:lx2,1:lx3,3))
    allocate(self%etheta, self%ephi, self%e1, self%e2, self%e3, mold=self%er)

    allocate(self%I(1:lx2,1:lx3))
    allocate(self%Bmag(1:lx1,1:lx2,1:lx3))
!    allocate(self%g1, self%g2, self%g3, self%r, self%theta, self%phi, self%alt, self%glon, self%glat, mold=self%Bmag)
    allocate(self%g1, self%g2, self%g3, mold=self%Bmag)

    ! coordinate-related quantities will retain ghost cells
    allocate(self%r(-1:lx1+2,-1:lx2+2,-1:lx3+2))
    allocate(self%theta, self%phi, self%alt, self%glon, self%glat, mold=self%r)

    !! by default these are not use so we allocate on request in set_geographici()
    !! some location edges are required; only for workers tho
    !allocate(self%gloni(1:lx1+1,1:lx2+1,1:lx3+1))
    !allocate(self%glati,self%alti, mold=self%gloni)

    self%coord_alloc_status=.true.
  end subroutine init_storage


  !> allocate space for root-only grid quantities
  subroutine init_storage_root(self)
    class(curvmesh), intent(inout) :: self
    integer :: lx1,lx2all,lx3all

    !fixme:  allocate root storage here; maybe check myid==0???  Need some protection so this does not get
    !         called from a worker...  Assume for now the calling procedure provides that.

    lx1=self%lx1
    lx2all=self%lx2all
    lx3all=self%lx3all

    if(lx1 < 1) error stop 'meshobj:init_storage_root: lx1 must be strictly positive'
    if(lx2all < 1) error stop 'meshobj:init_storage_root: lx2all must be strictly positive'
    if(lx3all < 1) error stop 'meshobj:init_storage_root: lx3all must be strictly positive'

    if (self%coord_alloc_status_root) error stop 'attempting to allocate space for root-only arrays when they already exist!'
    allocate(self%h1all(-1:lx1+2,-1:lx2all+2,-1:lx3all+2))
    allocate(self%h2all, self%h3all, mold=self%h1all)

    allocate(self%h1x1iall(1:lx1+1,1:lx2all,1:lx3all))
    allocate(self%h2x1iall, self%h3x1iall, mold=self%h1x1iall)

    allocate(self%h1x2iall(1:lx1,1:lx2all+1,1:lx3all))
    allocate(self%h2x2iall, self%h3x2iall, mold=self%h1x2iall)

    allocate(self%h1x3iall(1:lx1,1:lx2all,1:lx3all+1))
    allocate(self%h2x3iall, self%h3x3iall, mold=self%h1x3iall)

!    allocate(self%rall(1:lx1,1:lx2all,1:lx3all))
!    allocate(self%thetaall, self%phiall, self%altall, self%Bmagall, self%glonall, mold=self%rall)

    allocate(self%rall(-1:lx1+2,-1:lx2all+2,-1:lx3all+2))
    allocate(self%thetaall, self%phiall, self%altall, self%glonall, mold=self%rall)

    allocate(self%Bmagall(1:lx1,1:lx2all,1:lx3all))

    self%coord_alloc_status_root=.true.
  end subroutine init_storage_root


  !> this sets the geographic center of the mesh; must be called prior to make()
  subroutine set_center(self,glon,glat)
    class(curvmesh), intent(inout) :: self
    real(wp), intent(in) :: glon,glat

    self%glonctr=glon
    self%glatctr=glat
  end subroutine set_center


  !> determine what type of grid we have
  subroutine calc_gridflag(self)
    class(curvmesh), intent(inout) :: self

    integer :: lx1,lx2,lx3

    ! error checking
    if (.not. self%geog_set_status) error stop ' attempting to compute gridflag prior to setting geographic coordinates!'

    ! sizes
    lx1=self%lx1
    lx2=self%lx2
    lx3=self%lx3

    ! grid type, assumes that each worker x1 arrays span the entire altitude range of the grid
    if (abs(self%alt(1,1,1)-self%alt(lx1,1,1)) < 100e3_wp) then    !closed dipole grid
      self%gridflag=0
    else if (self%alt(1,1,1) > self%alt(2,1,1)) then            !open dipole grid with inverted structure wrt altitude
      self%gridflag=1
    else                                                      !something different (viz. non-inverted - lowest altitudes at the logical bottom of the grid)
      self%gridflag=2
    end if
  end subroutine calc_gridflag


  !> compute the number of null grid points and their indices for later use
  pure subroutine calc_inull(self)
    class(curvmesh), intent(inout) :: self
    integer :: lx1,lx2,lx3
    integer :: icount,ix1,ix2,ix3
    real(wp), dimension(:,:,:), allocatable :: alttmp

    ! error checking, we require the geographic coords. before this is done
    if (.not. self%geog_set_status) error stop ' attempting to compute null points prior to setting geographic coordinates!'

    ! sizes
    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    ! set null points for this simulation
    allocate(self%nullpts(1:lx1,1:lx2,1:lx3))
    self%nullpts=.false.
    allocate(alttmp(lx1,lx2,lx3))
    alttmp=self%alt(1:lx1,1:lx2,1:lx3)
    where (alttmp < 50e3_wp)
      self%nullpts=.true.
    end where
    deallocate(alttmp)

    ! count needed storage for null indices
    self%lnull=0;
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          if(self%nullpts(ix1,ix2,ix3)) self%lnull=self%lnull+1
        end do
      end do
    end do
    allocate(self%inull(1:self%lnull,1:3))

    ! store null indices
    icount=1
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          if(self%nullpts(ix1,ix2,ix3)) then
            self%inull(icount,:)=[ix1,ix2,ix3]
            icount=icount+1
          end if
        end do
      end do
    end do

    self%null_alloc_status=.true.
  end subroutine calc_inull


  !> compute geographic coordinates of all grid points
  pure subroutine calc_geographic(self)
    class(curvmesh), intent(inout) :: self

    ! FIXME: error checking?
    ! do the geographic conversion
    call geomag2geog(self%phi,self%theta,self%glon,self%glat)
    self%alt=r2alt(self%r)
    self%geog_set_status=.true.
  end subroutine calc_geographic


  !> this isn't used by computational parts of gemini but may be need for certain types of 3D visualization
  !    as such we only allocate and calculate interface geographic locations on demand and not be default
  subroutine calc_geographici(self)
    class(curvmesh), intent(inout) :: self
    real(wp), dimension(:,:,:), allocatable :: ri,thetai,phispheri

    if (.not. self%geogi_set_status) then 
      allocate(self%alti(1:self%lx1+1,1:self%lx2+1,1:self%lx3+1))
      allocate(self%gloni,self%glati, mold=self%alti) 

      allocate(ri,thetai,phispheri, mold=self%alti)

      call self%native2ECEFspher(self%glonctr,self%glatctr,self%x1i,self%x2i,self%x3i,ri,thetai,phispheri)
      !call geomag2geog(phispheri,thetai,self%gloni,self%glati) 
      !self%alti=r2alt(ri)

      self%alti(:,:,:)=ri(:,:,:)*cos(thetai(:,:,:))    ! z
      self%gloni(:,:,:)=ri(:,:,:)*sin(thetai(:,:,:))*cos(phispheri(:,:,:))    ! x
      self%glati(:,:,:)=ri(:,:,:)*sin(thetai(:,:,:))*sin(phispheri(:,:,:))    ! y

      deallocate(ri,thetai,phispheri)
      self%geogi_set_status=.true.
    else
      !print*, 'WARNING:  geographic locations of cell edges already calculated!'
      return
    end if
  end subroutine calc_geographici


  !> procedure to compute (but not store - external arrays provided as input) and geographic coordinate unit vectors
  !    This works on a full spatial arrays worth of data.
  subroutine calc_unitvec_geo(self,ealt,eglon,eglat)
    class(curvmesh), intent(in) :: self
    real(wp), dimension(:,:,:,:), intent(out) :: ealt,eglon,eglat
    integer :: lx1,lx2,lx3,ix1,ix2,ix3
    real(wp) :: thetagg,phigg    ! geographic spherical coords
    real(wp), dimension(3,3) :: Rgg2gm
    real(wp), dimension(3,1) :: ehere,ehererot

    if ( .not. self%geog_set_status) then
      error stop 'geographic coords. must be set prior to computing unit vectors'
    endif
    ! sizes
    lx1=self%lx1; lx2=self%lx2; lx3=self%lx3

    ! rotate geographic Cartesian ECEF into geomagnetic Cartesian ECEF
    Rgg2gm=rotgg2gm()
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          ! spherical geographic coords
          thetagg=pi/2-self%glat(ix1,ix2,ix3)*pi/180
          phigg=self%glon(ix1,ix2,ix3)*pi/180

          ! altitude unit vector
          ehere(1:3,1)=er_spherical(thetagg,phigg)
          ehererot=matmul(Rgg2gm,ehere)
          ealt(ix1,ix2,ix3,1:3)=ehererot(1:3,1)

          ! latitude unit vector
          ehere(1:3,1) = -1*etheta_spherical(thetagg,phigg)
          ehererot=matmul(Rgg2gm,ehere)
          eglat(ix1,ix2,ix3,1:3)=ehererot(1:3,1)

          ! longitude
          ehere(1:3,1) = ephi_spherical(thetagg,phigg)
          ehererot=matmul(Rgg2gm,ehere)
          eglon(ix1,ix2,ix3,1:3)=ehererot(1:3,1)
        end do
      end do
    end do
  end subroutine calc_unitvec_geo


  !> write the size of the grid to a file
  subroutine writesize(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation

    IDstr=' '
    write(IDstr,'(i1)') ID
    call hf%open(path//'/simsize'//IDstr//'.h5', action='w')
    call hf%write('/lx1',self%lx1)
    call hf%write('/lx2',self%lx2)
    call hf%write('/lx3',self%lx3)
    call hf%close()
  end subroutine writesize


  !> write grid coordinates (curvilinear only) to a file
  subroutine writegrid(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation

    ! at a minimum we must have allcated the coordinate arrays
    if (.not. self%xi_alloc_status) error stop ' attempting to write coordinate arrays when they have not been set yet!'

    IDstr=' '
    write(IDstr,'(i1)') ID
    call self%writesize(path,ID)
    call hf%open(path//'/simgrid'//IDstr//'.h5', action='w')
    call hf%write('/x1',self%x1)
    call hf%write('/x2',self%x2)
    call hf%write('/x3',self%x3)
    call hf%close()
  end subroutine writegrid


  subroutine writegridall(self,path,ID)
    class(curvmesh), intent(in) :: self
    character(*), intent(in) :: path
    integer, intent(in) :: ID
    type(hdf5_file) :: hf
    character(:), allocatable :: IDstr    ! use auto-allocation
    real(wp), dimension(:,:,:), allocatable :: realnullpts

    ! at a minimum we must have allcated the coordinate arrays
    if (.not. self%xi_alloc_status) error stop ' attempting to write coordinate arrays when they have not been set yet!'

    IDstr=' '
    write(IDstr,'(i1)') ID
    call self%writesize(path,ID)
    call hf%open(path//'/simgrid'//IDstr//'.h5', action='w')
    call hf%write('/x1',self%x1)
    call hf%write('/x1i',self%x1i)
    call hf%write('/dx1b',self%dx1)
    call hf%write('/dx1h',self%dx1i)

    call hf%write('/x2',self%x2)
    call hf%write('/x2i',self%x2i)
    call hf%write('/dx2b',self%dx2)
    call hf%write('/dx2h',self%dx2i)

    call hf%write('/x3',self%x3)
    call hf%write('/x3i',self%x3i)
    call hf%write('/dx3b',self%dx3)
    call hf%write('/dx3h',self%dx3i)

    call hf%write('/h1',self%h1)
    call hf%write('/h2',self%h2)
    call hf%write('/h3',self%h3)

    call hf%write('/h1x1i',self%h1x1i)
    call hf%write('/h2x1i',self%h2x1i)
    call hf%write('/h3x1i',self%h3x1i)

    call hf%write('/h1x2i',self%h1x2i)
    call hf%write('/h2x2i',self%h2x2i)
    call hf%write('/h3x2i',self%h3x2i)

    call hf%write('/h1x3i',self%h1x3i)
    call hf%write('/h2x3i',self%h2x3i)
    call hf%write('/h3x3i',self%h3x3i)

    call hf%write('/gx1',self%g1)
    call hf%write('/gx2',self%g2)
    call hf%write('/gx3',self%g3)

    call hf%write('/alt',self%alt)
    call hf%write('/glon',self%glon)
    call hf%write('/glat',self%glat)

    call hf%write('/Bmag',self%Bmag)
    call hf%write('/I',self%I)
    allocate(realnullpts(1:self%lx1,1:self%lx2,1:self%lx3))
    realnullpts=0.0
    where (self%nullpts)
      realnullpts = 1
    end where
    call hf%write('/nullpts',realnullpts)
    deallocate(realnullpts)

    call hf%write('/e1',self%e1)
    call hf%write('/e2',self%e2)
    call hf%write('/e3',self%e3)

    call hf%write('/er',self%er)
    call hf%write('/etheta',self%etheta)
    call hf%write('/ephi',self%ephi)

    call hf%write('/r',self%r)
    call hf%write('/theta',self%theta)
    call hf%write('/phi',self%phi)

    call hf%close()
  end subroutine writegridall


  !> deallocate space associated with pointers and set appropriate flags
  subroutine dissociate_pointers(self)
    class(curvmesh), intent(inout) :: self

    ! deallocation statements here; always check allocation status flags first...
    if (self%xi_alloc_status) then
      deallocate(self%x1,self%x2,self%x3,self%x2all,self%x3all)    ! these are from set_coords
      self%xi_alloc_status=.false.
    end if
    if (self%dxi_alloc_status) then                                  ! from calc_coord_diffs
      deallocate(self%dx1,self%x1i,self%dx1i)
      deallocate(self%dx2,self%x2i,self%dx2i)
      deallocate(self%dx3,self%x3i,self%dx3i)
      self%dxi_alloc_status=.false.
    end if
    if (self%dxi_alloc_status_root) then
      deallocate(self%dx2all,self%x2iall,self%dx2iall)
      deallocate(self%dx3all,self%x3iall,self%dx3iall)
      self%dxi_alloc_status_root=.false.
    end if
    if (self%difflen_alloc_status) then
      deallocate(self%dl1i,self%dl2i,self%dl3i)    ! from calc_difflengths
      self%difflen_alloc_status=.false.
    end if

    ! coordinate-specific arrays set by type extensions
    if (self%coord_alloc_status) then
      deallocate(self%er,self%etheta,self%ephi)
      deallocate(self%e1,self%e2,self%e3)
      deallocate(self%r,self%theta,self%phi)
      deallocate(self%h1,self%h2,self%h3)
      deallocate(self%h1x1i,self%h2x1i,self%h3x1i)
      deallocate(self%h1x2i,self%h2x2i,self%h3x2i)
      deallocate(self%h1x3i,self%h2x3i,self%h3x3i)
      deallocate(self%g1,self%g2,self%g3)
      deallocate(self%Bmag,self%I)
      deallocate(self%alt,self%glon,self%glat)
      self%coord_alloc_status=.false.
      self%geog_set_status=.false.
    end if
    if (self%geogi_set_status) then
      deallocate(self%alti,self%gloni,self%glati)
      self%geogi_set_status=.false.
    end if
    if (self%coord_alloc_status_root) then
      deallocate(self%h1all,self%h2all,self%h3all)
      deallocate(self%h1x1iall,self%h2x1iall,self%h3x1iall)
      deallocate(self%h1x2iall,self%h2x2iall,self%h3x2iall)
      deallocate(self%h1x3iall,self%h2x3iall,self%h3x3iall)
      deallocate(self%rall,self%thetaall,self%phiall)
      deallocate(self%altall,self%Bmagall)
      deallocate(self%glonall)
      self%coord_alloc_status_root=.false.
      self%coord_set_status_root=.false.
    end if
    if (self%null_alloc_status) then
      deallocate(self%nullpts,self%inull)
      self%null_alloc_status=.false.
    end if
  end subroutine dissociate_pointers
  ! FIXME:  it may make sense to have a procedure to clear unit vectors here to conserve some memory
  ! FIXME:  also have some way to clear the fullgrid metric factors once they are dealt with?
end module meshobj
