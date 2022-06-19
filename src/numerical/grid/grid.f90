!> non-mpi related grid quantities are stored here for use by other modules.  This could arguably be an 
!    object...
module grid

use meshobj, only: curvmesh
use phys_consts, only: wp
use reader, only: get_simsize3

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occasionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
real(wp), protected :: glonctr=-720._wp,glatctr=-180._wp

private 
public :: lx1,lx2,lx3,lx2all,lx3all,gridflag, &
             get_grid3_coords_raw, get_grid3_coords_hdf5, get_grid3_coords_nc4, &
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,grid_size, &
             grid_from_extents, generate_worker_grid, ungenerate_worker_grid, &
             get_grid3_coords,read_size_gridcenter,detect_gridtype,set_size_gridcenter

interface ! readgrid_*.f90
  module subroutine get_grid3_coords_raw(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp), intent(out) :: glonctr,glatctr
  end subroutine get_grid3_coords_raw
  module subroutine get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp), intent(out) :: glonctr,glatctr
  end subroutine get_grid3_coords_hdf5
  module subroutine get_grid3_coords_nc4(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp), intent(out) :: glonctr,glatctr
  end subroutine get_grid3_coords_nc4
end interface

contains
  !> detect the type of grid that we are dealing with based solely on native coordinate values
  function detect_gridtype(x1,x2,x3) result(xtype)
    real(wp), dimension(-1:), intent(in) :: x1,x2,x3
    integer :: xtype
   
    if (maxval(abs(x2))<100) then
      print*, ' Detected dipole grid...'
      xtype=2
    else
      print*, 'Detected Cartesian grid...'
      xtype=1
    end if
  end function detect_gridtype


  !> Force a size and grid center location into module variables
  subroutine set_size_gridcenter(lx1in,lx2allin,lx3allin,glonctrin,glatctrin)
    integer, intent(in) :: lx1in,lx2allin,lx3allin
    real(wp), intent(in) :: glonctrin,glatctrin

    lx1=lx1in; lx2all=lx2allin; lx3all=lx3allin;
    glonctr=glonctrin; glatctr=glatctrin;
  end subroutine set_size_gridcenter


  !> Query the coordinates file and pull out the center geographic location for the entire grid (used for
  !    generation of Cartesian meshes and put in in a module-scope variable.
  subroutine read_size_gridcenter(indatsize,outdir)
    character(*), intent(in) :: indatsize,outdir
    real(wp), dimension(:), allocatable :: x1,x2all,x3all
    
    call get_simsize3(indatsize,lx1,lx2all,lx3all)
    allocate(x1(-1:lx2all+2),x2all(-1:lx2all+2),x3all(-1:lx3all+2))
    call get_grid3_coords(outdir,x1,x2all,x3all,glonctr,glatctr)
    deallocate(x1,x2all,x3all)
  end subroutine read_size_gridcenter


  !> Generate grid from a set of extents and sizes - e.g. similar to what is used in forestcalw.  input
  !    sizes should include ghost cells.  WARNING: this function will always just assume you are using a 
  !    local grid, i.e. one that doesn't need knowledge of the full grid extents!
  subroutine grid_from_extents(x1lims,x2lims,x3lims,lx1wg,lx2wg,lx3wg,x)
    real(wp), dimension(2), intent(in) :: x1lims,x2lims,x3lims
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    class(curvmesh), intent(inout) :: x
    integer :: ix1,ix2,ix3
    real(wp), dimension(:), allocatable :: x1,x2,x3

    ! error checking
    if (glatctr<-90._wp .or. glatctr>90._wp) then
      error stop ' grid_from_extents:  prior to calling must use read_size_gridcenter or set_size_gridcenter to assign &
                   module variables glonctr,glatctr'
    end if

    ! create temp space
    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays
    x1=[(x1lims(1) + (x1lims(2)-x1lims(1))/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lims(1) + (x2lims(2)-x2lims(1))/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lims(1) + (x3lims(2)-x3lims(1))/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)] 

    ! generate a subgrid from these
    call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)

    ! get rid of temp. arrays
    deallocate(x1,x2,x3)
  end subroutine grid_from_extents


  !> Generate a "worker" grid based soley on coordinate arrays, polymorphic grid object must already
  !    exist, i.e. already be allocated with some dynamic type.  Note that you can set x2all=x2 and 
  !    (or) x3all=x3 if you are only doing "local" grid operations in your GEMINI application, e.g. as with 
  !    trees-GEMINI.
  subroutine generate_worker_grid(x1,x2,x3,x2all,x3all,glonctr,glatctr,x)
    real(wp), dimension(:), intent(in) :: x1,x2,x3,x2all,x3all
    real(wp), intent(in) :: glonctr,glatctr 
    class(curvmesh), intent(inout) :: x

    ! Create the grid object
    call x%set_center(glonctr,glatctr)
    call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
    call x%init()                              ! allocate space for subgrid variables
    call x%make()                              ! fill auxiliary arrays
  
    call set_gridflag(x%gridflag)
  
    !> Set gravitational fields for module scope vars., use pointers to avoid duplicating data
    !call bind_grav_ptrs(x%g1,x%g2,x%g3)  
  end subroutine generate_worker_grid
  
  
  !> Force deallocation of grid data at least to the point where it can be "remade", e.g. for AMR-like operations
  subroutine ungenerate_worker_grid(x)
    class(curvmesh), intent(inout) :: x
  
    call x%dissociate_pointers()
  end subroutine ungenerate_worker_grid


  !> Read in native coordinates from a grid file
  subroutine get_grid3_coords(path,x1,x2all,x3all,glonctr,glatctr)
    character(*), intent(in) :: path
    real(wp), dimension(:), intent(inout) :: x1,x2all,x3all
    real(wp) :: glonctr,glatctr
  
    character(:), allocatable :: fmt
  
    fmt = path(index(path, '.', back=.true.) : len(path))
    select case (fmt)
      case ('.dat')
        call get_grid3_coords_raw(path,x1,x2all,x3all,glonctr,glatctr)
      case ('.h5')
        call get_grid3_coords_hdf5(path,x1,x2all,x3all,glonctr,glatctr)
      case ('.nc')
        call get_grid3_coords_nc4(path,x1,x2all,x3all,glonctr,glatctr)
      case default
        error stop 'grid:read:get_grid3: unknown grid format: ' // fmt
    end select
  
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
