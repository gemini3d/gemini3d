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

!! FIXME:  these must be stored in a patch user_data object!!!!
!real(wp), dimension(:,:,:), pointer, protected :: g1,g2,g3
!! so other modules can simply access gravitational field data

public :: lx1,lx2,lx3,lx2all,lx3all,gridflag, &
             get_grid3_coords_raw, get_grid3_coords_hdf5, get_grid3_coords_nc4, &
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,grid_size

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

contains    !! all we have are setter procedures + whatever mpi-independent stuff is in submodules
  !> Generate grid from a set of extents and sizes - e.g. similar to what is used in forestcalw.  input
  !    sizes should include ghost cells.  WARNING: this function will always just assume you are using a 
  !    local grid, i.e. one that doesn't need knowledge of the full grid extents!
  subroutine grid_from_extents(x1lower,x1upper,x2lower,x2upper,x3lower,x3upper,lx1wg,lx2wg,lx3wg,glonctr,glatctr,x)
    real(wp), intent(in) :: x1lower,x1upper,x2lower,x2upper,x3lower,x3upper
    integer, intent(in) :: lx1wg,lx2wg,lx3wg
    real(wp), intent(in) :: glonctr,glatctr
    class(curvmesh), intent(inout) :: x
    integer :: ix1,ix2,ix3
    real(wp), dimension(:), allocatable :: x1,x2,x3

    ! create temp space
    allocate(x1(lx1wg),x2(lx2wg),x3(lx3wg))

    ! make uniformly spaced coordinate arrays
    x1=[(x1lower + (x1upper-x1lower)/(lx1wg-1)*(ix1-1),ix1=1,lx1wg)]
    x2=[(x2lower + (x2upper-x2lower)/(lx2wg-1)*(ix2-1),ix2=1,lx2wg)]
    x3=[(x3lower + (x3upper-x3lower)/(lx3wg-1)*(ix3-1),ix3=1,lx3wg)] 

    ! generate a subgrid from these
    call generate_worker_grid(x1,x2,x3,x2,x3,glonctr,glatctr,x)

    ! release memory
    deallocate(x1,x2,x3)
  end subroutine grid_from_extents


  !> Generate a "worker" grid based soley on coordinate arrays, polymorphic grid object must already
  !    exist, i.e. already be allocated.  Note that you can set x2all=x2 and x3all=x3 if you are only
  !    doing "local" grid operations in your GEMINI application, e.g. as with trees-GEMINI.
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
