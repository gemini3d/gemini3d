!> non-mpi related grid quantities are stored here for use by other modules.  This could arguably be an
!    object...
module grid

use, intrinsic :: iso_c_binding, only: C_PTR,C_INT,c_loc,c_f_pointer,C_NULL_PTR
use meshobj, only: curvmesh
use meshobj_dipole, only: dipolemesh
use meshobj_cart, only: cartmesh
use phys_consts, only: wp
use reader, only: get_simsize3

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occasionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
real(wp), protected :: glonctr=-720._wp,glatctr=-180._wp
real(wp), dimension(2), protected :: x1lims,x2alllims,x3alllims
real(wp), dimension(:), pointer, protected :: x1=>null()
logical, protected :: flaglims=.false.
!!^ These variables will be shared between all workers/patches so they can be module-scope variables.

private
public :: lx1,lx2,lx3,lx2all,lx3all,gridflag,x1, &
             get_grid3_coords_hdf5, &
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,grid_size, &
             grid_from_extents, grid_internaldata_alloc, grid_internaldata_generate, &
             grid_internaldata_ungenerate, &
             get_grid3_coords,read_size_gridcenter,detect_gridtype,set_size_gridcenter, &
             meshobj_alloc, get_gridcenter, meshobj_dealloc, set_fullgrid_lims, &
             x1lims,x2alllims,x3alllims,get_x1coords,get_fullgrid_lims, alloc_x1coords, &
             read_grid, grid_check, grid_drift, calc_subgrid_size

             !, generate_worker_grid

interface ! readgrid_*.f90
  module subroutine get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp), intent(out) :: glonctr,glatctr
  end subroutine
end interface

interface !< grid_mpi.f90
module subroutine read_grid(indatsize,indatgrid,flagperiodic, x, xtype, xC)
  !1 read in grid and set subgrid sizes; total size must already be set in the grid module via grid_size().
  !!    this is only to be used when GEMINI is run using functionality that depends on fullgrid data like
  !!    potential solutions etc.
  character(*), intent(in) :: indatsize,indatgrid
  integer, intent(in) :: flagperiodic
  class(curvmesh), pointer, intent(inout) :: x
  integer(C_INT), intent(inout), optional :: xtype
  type(C_PTR), intent(inout), optional :: xC
end subroutine

module subroutine calc_subgrid_size(lx2all, lx3all)
  !! worker subgrid sizes; requires knowledge of mpi, though not any direct mpi calls
  integer, intent(in) :: lx2all, lx3all
end subroutine

module subroutine grid_drift(x,E02,E03,v2grid,v3grid)
  !1 Compute grid drift speed; requires that we exchange some data through mpi
  !! Compute the speed the grid is moving at given a background electric field
  class(curvmesh), intent(in) :: x
  reaL(wp), dimension(:,:,:), intent(in) :: E02,E03
  real(wp), intent(inout) :: v2grid,v3grid
  !! intent(out)
end subroutine
end interface

interface !< check.f90
  module subroutine grid_check(x)
    class(curvmesh), intent(in) :: x
  end subroutine
end interface


!! some overloading for situations needing vs. not needing an allocation step
!interface grid_from_extents
!  module procedure grid_from_extents_noalloc,grid_from_extents_alloc
!end interface

contains
  !> detect the type of grid that we are dealing with based solely on native coordinate values
  function detect_gridtype(x1,x2,x3) result(xtype)
    real(wp), dimension(-1:), intent(in) :: x1,x2,x3
    integer :: xtype

    if (maxval(abs(x2))<1000) then
      print '(a)', 'Detected dipole grid...'
      xtype=2
    else
      print '(a)', 'Detected Cartesian grid...'
      xtype=1
    end if
  end function detect_gridtype


  !> Force a size and grid center location into module variables, if desired.  In general some other method
  !    should be used like read_size_gridcenter().
  subroutine set_size_gridcenter(lx1in,lx2allin,lx3allin,glonctrin,glatctrin)
    integer, intent(in) :: lx1in,lx2allin,lx3allin
    real(wp), intent(in) :: glonctrin,glatctrin

    lx1=lx1in; lx2all=lx2allin; lx3all=lx3allin;
    glonctr=glonctrin; glatctr=glatctrin;
  end subroutine set_size_gridcenter


  !> retrieve the grid center location from module
  subroutine get_gridcenter(glonctrout,glatctrout)
    real(wp), intent(inout) :: glonctrout,glatctrout

    glonctrout=glonctr
    glatctrout=glatctr
  end subroutine get_gridcenter


  subroutine alloc_x1coords(lx1)
    integer, intent(in) :: lx1

    allocate(x1(-1:lx1+2))
  end subroutine alloc_x1coords


  !> retrieve the x1 coordinate array from module
  subroutine get_x1coords(x1out)
    real(wp), dimension(:), pointer, intent(inout) :: x1out

    if (associated(x1)) then
      x1out=>x1
    else
      print*, 'WARNING - grid:get_x1coords - x1 variable is not allocated/set, check call sequence.'
    end if
  end subroutine get_x1coords


  !> Query the coordinates file and pull out the center geographic location for the entire grid (used for
  !    generation of Cartesian meshes and put in in a module-scope variable.
  subroutine read_size_gridcenter(indatsize,indatgrid)
    character(*), intent(in) :: indatsize,indatgrid
    real(wp), dimension(:), allocatable :: x2all,x3all

    call get_simsize3(indatsize,lx1,lx2all,lx3all)
    call alloc_x1coords(lx1)    ! module-scope
    allocate(x2all(-1:lx2all+2),x3all(-1:lx3all+2))
    call get_grid3_coords(indatgrid,x1,x2all,x3all,glonctr,glatctr)
    ! FIXME: should store min/max here; can be used to detect whether we are on the global boundary.  We'd also need
    !   to add this data from other grid creation interfaces in the grid_mpi.f90 module.
    call set_fullgrid_lims(x2all,x3all)
    deallocate(x2all,x3all)
  end subroutine read_size_gridcenter


  !> set the fullgrid limit variables in the model, e.g. for detecting if we are on a global boundary
  subroutine set_fullgrid_lims(x2all,x3all)
    real(wp), dimension(-1:) :: x2all
    real(wp), dimension(-1:) :: x3all

    x1lims=[x1(1),x1(lx1)]
    x2alllims=[x2all(1),x2all(lx2all)]
    x3alllims=[x3all(1),x3all(lx3all)]
    flaglims=.true.
  end subroutine set_fullgrid_lims


  !> return the extents of the FULL GRID
  subroutine get_fullgrid_lims(x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax)
    real(wp), intent(inout) :: x1min,x1max,x2allmin,x2allmax,x3allmin,x3allmax

    if (flaglims) then
      x1min=x1lims(1); x1max=x1lims(2);
      x2allmin=x2alllims(1); x2allmax=x2alllims(2);
      x3allmin=x3alllims(1); x3allmax=x3alllims(2);
    else
      error stop 'grid:get_fullgrid_lims - attempt to retrieve grid limits when not set...'
    end if
  end subroutine get_fullgrid_lims


  !! FIXME: deprecated?
  !> Generate grid from a set of extents and sizes - e.g. similar to what is used in forestcalw.  input
  !    sizes should include ghost cells.  WARNING: this function will always just assume you are using a
  !    local grid, i.e. one that doesn't need knowledge of the full grid extents!  This requires that
  !    the grid type/class already be defined.
  subroutine grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x
    integer :: ix1,ix2,ix3
    real(wp), dimension(:), allocatable :: x1,x2,x3

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop 'ERROR:grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter' // &
        'to assign module variables glonctr,glatctr'
    end if

    ! create temp space
    !allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))

    ! see if we need to allocate x1 module variable or if it is already set up
    if (.not. allocated(x1)) then
      allocate(x1(lx1wg))
      x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]    ! usually x1 is nonuniform and already set up here but if not allocate and set using uniform forcing.  If already allocated then the limits will be ignored.
    end if
    ! make uniformly spaced coordinate arrays
    allocate(x2(lx2wg),x3(lx3wg))
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]

    !call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)
    call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)
    call grid_internaldata_generate(x)

    ! get rid of temp. arrays
    deallocate(x2,x3)
  end subroutine grid_from_extents


  !> use a ghostgrid object to retrieve data about ghost cell geographic locations for use, e.g., in generating
  !    vtu output files for ForestGEMINI.  The output arrays first two dimensions represent the variations with
  !    location on the boundary whereas the 3rd dim. is the component of position vector ordered at glon,glat,alt
!  subroutine ghost_location_generate(x,gcoordsx1max,gcoordsx2max,gcoordsx3max)
!    class(curvmesh), intent(in) :: x
!    real(wp), dimension(:,:,:), intent(inout) :: gcoordsx1max
!    real(wp), dimension(:,:,:), intent(inout) :: gcoordsx2max
!    real(wp), dimension(:,:,:), intent(inout) :: gcoordsx3max
!    class(curvmesh), allocatable :: xghost
!    real(wp) :: x1max,x2max,x3max
!
!    ! error check for size
!    if (.not. (size(gcoordsx1max,1)==x%lx2 .and. size(gcoordsx1max,2)==x%lx3 .and. size(gcoordsx1max,3)==3) ) then
!      error stop 'ghost_location_generate:  bad array size in gcoordsx1max spec'
!    end if
!    if (.not. (size(gcoordsx2max,1)==x%lx1 .and. size(gcoordsx2max,2)==x%lx3 .and. size(gcoordsx2max,3)==3) ) then
!      error stop 'ghost_location_generate:  bad array size in gcoordsx2max spec'
!    end if
!    if (.not. (size(gcoordsx3max,1)==x%lx1 .and. size(gcoordsx3max,2)==x%lx2 .and. size(gcoordsx3max,3)==3) ) then
!      error stop 'ghost_location_generate:  bad array size in gcoordsx1max spec'
!    end if
!
!    ! create a "ghost grid" object that has internal cells corresponding to ghost locations in actual grid
!    select type (x)
!      type is (dipolemesh)
!        allocate(dipolemesh::xghost)
!      type is (cartmesh)
!        allocate(cartmesh::xghost)
!      default
!        error stop 'ghost_location_generate:  could not determine type of source mesh'
!    end select
!
!    ! create coordinate for the ghostgrid just corresponding to the 3 slices of ghostcells we need
!    allocate(x1tmp(5),x2tmp(5),x3tmp(5))
!    x1tmp=[x%x1(lx1-1),x%x1(lx1),x%x1(lx1+1),x%x1(lx1+2),x%x1(lx1+2)+0.1]
!    call grid_internaldata_alloc(x1tmp,x%x2,x%x3,x%x2,x%x3,0.0,0.0,xghost)
!    call grid_internaldata_generate(xghost)
!    gcoordsx1max(1:lx2,1:lx3,1)=xghost%glon(1,1:lx2,1:lx3)
!    gcoordsx1max(1:lx2,1:lx3,2)=xghost%glat(1,1:lx2,1:lx3)
!    gcoordsx1max(1:lx2,1:lx3,3)=xghost%alt(1,1:lx2,1:lx3)
!    call grid_internaldata_ungenerate(xghost)
!
!    x2tmp=[x%x2(lx2-1),x%x2(lx2),x%x2(lx2+1),x%x2(lx2+2),x%x2(lx2+2)+0.1]
!    call grid_internaldata_alloc(x%x1,x2tmp,x%x3,x2tmp,x%x3,0.0,0.0,xghost)
!    call grid_internaldata_generate(xghost)
!    gcoordsx1max(1:lx1,1:lx3,1)=xghost%glon(1:lx1,1,1:lx3)
!    gcoordsx1max(1:lx1,1:lx3,2)=xghost%glat(1:lx1,1,1:lx3)
!    gcoordsx1max(1:lx1,1:lx3,3)=xghost%alt(1:lx1,1,1:lx3)
!    call grid_internaldata_ungenerate(xghost)
!
!    x3tmp=[x%x3(lx3-1),x%x3(lx3),x%x3(lx3+1),x%x3(lx3+2),x%x3(lx3+2)+0.1]
!    call grid_internaldata_alloc(x%x1,x%x2,x3tmp,x%x2,x3tmp,0.0,0.0,xghost)
!    call grid_internaldata_generate(xghost)
!    gcoordsx1max(1:lx1,1:lx2,1)=xghost%glon(1:lx1,1:lx2,1)
!    gcoordsx1max(1:lx1,1:lx2,2)=xghost%glat(1:lx1,1:lx2,1)
!    gcoordsx1max(1:lx1,1:lx2,3)=xghost%alt(1:lx1,1:lx2,1)
!    call grid_internaldata_ungenerate(xghost)
!
!    deallocate(x1tmp,x2tmp,x3tmp)
!    deallocate(xghost)
!  end subroutine ghost_location_generate


!  ! FIXME: split into grid_alloc and generate_worker_grid; add interfaces for both to libgemini and C
!  !> this version additionally allocates the input argument, which is now a pointer
!  subroutine grid_from_extents_alloc(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x,xtype,xC)
!    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
!    integer, intent(in) :: lx1wg,lx2wg,lx3wg
!    class(curvmesh), intent(inout), pointer :: x
!    integer, intent(inout) :: xtype
!    type(c_ptr), intent(inout) :: xC
!    integer :: ix1,ix2,ix3
!    real(wp), dimension(:), allocatable :: x1,x2,x3
!
!    ! error checking
!    if (glatctr<-90._wp .or. glatctr>90._wp) then
!      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
!                   module variables glonctr,glatctr'
!    end if
!
!    ! create temp space
!    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))
!
!    ! make uniformly spaced coordinate arrays
!    x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
!    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
!    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)]
!
!    ! generate a subgrid from these
!    !if (present(xtype) .and. present(xC)) then   ! optional arguments confuses overloading for some types of calls :/
!      call meshobj_alloc(x1,x2,x3,x,xtype,xC)
!    !else
!    !  call meshobj_alloc(x1,x2,x3,x2,x3,x)
!    !end if
!    call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)
!    ! call grid_internaldata_alloc(x1,x2,x3,x2,x3,glonctr,glatctr,x)
!    ! call grid_internaldata_generate(x)
!
!    ! get rid of temp. arrays
!    deallocate(x1,x2,x3)
!  end subroutine grid_from_extents_alloc


  !> Find the type of the grid and allocate the correct type/class, return a C pointer if requested
  !    via optional arguments
  subroutine meshobj_alloc(x1,x2,x3,x,xtype,xC)
    real(wp), dimension(:), intent(in) :: x1,x2,x3
    class(curvmesh), pointer, intent(inout) :: x
    integer(C_INT), intent(inout), optional :: xtype
    type(C_PTR), intent(inout), optional :: xC
    integer :: gridtype
    type(cartmesh), pointer :: xcart
    type(dipolemesh), pointer :: xdipole

    !! allocate and read correct grid type
    gridtype=detect_gridtype(x1,x2,x3)
    select case (gridtype)
      case (2)
        !allocate(dipolemesh::x)
        allocate(xdipole)
        x=>xdipole
        if (present(xC) .and. present(xtype)) then
          xC = c_loc(xdipole)
          xtype = gridtype
        end if
      case (1)
        !allocate(cartmesh::x)
        allocate(xcart)
        x=>xcart
        if (present(xC) .and. present(xtype)) then
          xC = c_loc(xcart)
          xtype = gridtype
        end if
      case default
        error stop 'grid:meshobj_alloc - Unable to identify grid type'
    end select
  end subroutine


  !> deallocate mesh class
  subroutine meshobj_dealloc(x,xtype,xC)
    class(curvmesh), pointer, intent(inout) :: x
    integer(C_INT), intent(inout), optional :: xtype
    type(C_PTR), intent(inout), optional :: xC


    deallocate(x)
    x=>null()
    xC=C_NULL_PTR
    !xtype=-1     ! an issue here is that the grid may be detected as "bad" even though deallocated correctly
  end subroutine meshobj_dealloc


!  !> Generate a "worker" grid based on coordinate arrays and grid center, polymorphic grid object must already
!  !    exist, i.e. already be allocated with some dynamic type.  Note that you can set x2all=x2 and
!  !    (or) x3all=x3 if you are only doing "local" grid operations in your GEMINI application, e.g. as with
!  !    trees-GEMINI.  The dynamic type of x must be set prior to calling this function; this can be
!  !    accomplished e.g. through a wrapper
!  subroutine generate_worker_grid(x1,x2,x3,x2all,x3all,glonctr,glatctr,x)
!    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
!    real(wp), intent(in) :: glonctr,glatctr
!    class(curvmesh), intent(inout) :: x
!
!    ! Create the grid object
!    call x%set_center(glonctr,glatctr)
!    call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
!    call x%init()                              ! allocate space for subgrid variables
!    call x%make()                              ! fill auxiliary arrays
!
!    call set_gridflag(x%gridflag)
!  end subroutine generate_worker_grid


  !> Trigger allocation of grid class internal data once the class itself has been allocated and typed
  subroutine grid_internaldata_alloc(x1,x2,x3,x2all,x3all,glonctr,glatctr,x)
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
    real(wp), intent(in) :: glonctr,glatctr
    class(curvmesh), intent(inout) :: x

    lx1=size(x1)-4
    lx2=size(x2)-4
    lx3=size(x3)-4
    call x%set_center(glonctr,glatctr)         ! set center location for grid (in case used, e.g. for Cartesian)
    call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
    call x%init()                              ! allocate space for subgrid variables
  end subroutine grid_internaldata_alloc


  !> Trigger a generation of all grid internal data
  subroutine grid_internaldata_generate(x)
    class(curvmesh), intent(inout) :: x

    call x%make()                              ! trigger generation of all internal data arrays
    call set_gridflag(x%gridflag)              ! set module variable to match the type stored in the grid class
  end subroutine


  !> Force deallocation of grid data at least to the point where it can be "remade", e.g. for AMR-like operations
  subroutine grid_internaldata_ungenerate(x)
    class(curvmesh), intent(inout) :: x

    call x%dissociate_pointers()
  end subroutine grid_internaldata_ungenerate


  !> Read in native coordinates from a grid file
  subroutine get_grid3_coords(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1
    real(wp), dimension(:), intent(inout) :: x2all,x3all
    real(wp) :: glonctr,glatctr

    call get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)

    if(size(x1) < 1) error stop 'grid:get_grid3_coords: size(x1) must be strictly positive'
    if(size(x2all) < 1) error stop 'grid:get_grid3_coords: size(x2all) must be strictly positive'
    if(size(x3all) < 1) error stop 'grid:get_grid3_coords: size(x3all) must be strictly positive'
  end subroutine get_grid3_coords


  subroutine set_total_grid_sizes(lx1in,lx2allin,lx3allin)
    integer, intent(in) :: lx1in,lx2allin,lx3allin

    lx1=lx1in; lx2all=lx2allin; lx3all=lx3allin;
  end subroutine set_total_grid_sizes


  subroutine set_subgrid_sizes(lx2in,lx3in)
    integer, intent(in) :: lx2in,lx3in

    lx2=lx2in; lx3=lx3in;
  end subroutine set_subgrid_sizes


  subroutine set_gridflag(gridflagin)
    integer, intent(in) :: gridflagin

    gridflag=gridflagin
  end subroutine set_gridflag


!  subroutine bind_grav_ptrs(g1in,g2in,g3in)
!    real(wp), dimension(:,:,:), pointer, intent(in) :: g1in,g2in,g3in
!
!    g1=>g1in; g2=>g2in; g3=>g3in
!  end subroutine bind_grav_ptrs


  subroutine grid_size(indatsize)
  !! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)
    character(*), intent(in) :: indatsize

    call get_simsize3(indatsize, lx1, lx2all, lx3all)
    print *, 'grid_size_root: full grid size:  ',lx1,lx2all,lx3all
    call set_total_grid_sizes(lx1,lx2all,lx3all)    !! set module global sizes for use on other contexts
  end subroutine grid_size
end module grid
