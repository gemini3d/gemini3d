!> non-mpi related grid quantities are stored here for use by other modules.  This could arguably be an 
!    object...
module grid

use phys_consts, only: wp
use reader, only: get_simsize3
use meshobj, only: curvmesh
use meshobj_dipole, only: dipolemesh
use meshobj_cart, only: cartmesh

integer, protected :: lx1,lx2,lx3,lx2all,lx3all
!! this is a useful shorthand for most program units using this module,
!! occasionally a program unit needs to define its own size in which case an only statement
!! is required when using this module.
integer, protected :: gridflag
!! for cataloguing the type of grid that we are using, open, closed, inverted, etc.  0 - closed dipole, 1 - inverted open, 2 - standard open.
real(wp), dimension(:,:,:), pointer, protected :: g1,g2,g3
!! so other modules can simply access gravitational field data

public :: lx1,lx2,lx3,lx2all,lx3all,gridflag, &
             get_grid3_coords_raw, get_grid3_coords_hdf5, get_grid3_coords_nc4, &
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,bind_grav_ptrs,grid_size

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

  subroutine bind_grav_ptrs(g1in,g2in,g3in)
    real(wp), dimension(:,:,:), pointer, intent(in) :: g1in,g2in,g3in

    g1=>g1in; g2=>g2in; g3=>g3in
  end subroutine bind_grav_ptrs

  subroutine grid_size(indatsize)
  !! CHECK THE SIZE OF THE GRID TO BE LOADED AND SET SIZES IN THIS MODULE (NOT IN STRUCTURE THOUGH)
    character(*), intent(in) :: indatsize
  
    call get_simsize3(indatsize, lx1, lx2all, lx3all)
    print *, 'grid_size_root: full grid size:  ',lx1,lx2all,lx3all
    call set_total_grid_sizes(lx1,lx2all,lx3all)    !! set module global sizes for use on other contexts
  end subroutine grid_size


  !> This performs a local grid setup given lower and upper bounds of the local grid and a number of grid points; this will assume
  !    uniform grid spacing
  subroutine setup_grid_uniform(bounds,flagperiodic,x)
    real(wp), dimension(6), intent(in) :: bounds
    integer, intent(in) :: flagperiodic
    class(curvmesh), allocatable, intent(inout) :: x
    real(wp) :: x1lower,x1upper,x2lower,x2upper,x3lower,x3upper
    real(wp), dimension(:), allocatable :: x1,x2,x3
    integer :: ix1,ix2,ix3
    real(wp) :: dx1,dx2,dx3  !grid step sizes

    ! convenience vars
    x1lower=bounds(1); x1upper=bounds(2);
    x2lower=bounds(3); x2upper=bounds(4);
    x3lower=bounds(5); x3upper=bounds(6);

    ! arrays to pass to grid construtor
    allocate(x1(-1:lx1+2),x2(-1:lx2+2),x3(-1:lx3+2))
    dx1=(x1upper-x1lower)/(lx1-1)
    do ix1=-1,lx1+2
      x1(ix1)=x1lower+(ix1-1)*dx1
    end do
    dx2=(x2upper-x2lower)/(lx2-1)
    do ix2=-1,lx2+2
      x2(ix2)=x2lower+(ix2-1)*dx2
    end do
    dx3=(x3upper-x3lower)/(lx3-1)
    do ix3=-1,lx3+2
      x3(ix3)=x3lower+(ix3-1)*dx3
    end do

    ! call the creation procedure
    call create_local_grid(x1,x2,x3,flagperiodic,x)
    deallocate(x1,x2,x3)
  end subroutine setup_grid_uniform


  !> Deallocate and (re)create the grid object based on input coordinate arrays.  Note that this will not force uniformity of geographic
  !    locations in x3 as needed for a constant background thermosphere.  
  subroutine create_local_grid(x1,x2,x3,flagperiodic,x)
    real(wp), dimension(-1:), intent(in) :: x1
    real(wp), dimension(-1:), intent(in) :: x2
    real(wp), dimension(-1:), intent(in) :: x3
    integer, intent(in) :: flagperiodic
    class(curvmesh), allocatable, intent(inout) :: x

    ! force existing grid structure to be deallocated and reallocated
    if (allocated(x)) then
      print*, ' Deallocating grid...'
      deallocate(x)
    end if
    
    ! right now we just have Cartesian and dipole so it's easy to detect based on x2
    if (maxval(abs(x2))<100) then
      print*, ' Allocating dipole grid...'
      allocate(dipolemesh::x)
    else
      print*, ' Allocating Cartesian grid...'
      allocate(cartmesh::x)
    end if

    ! generate grid, 
    ! FIXME: need to overload for subgrid v. root
    call x%set_local_coords(x1,x2,x3)    ! set primitive coordinates
    call x%init()                  ! allocate object arrays
    call x%make()                  ! compute metric factors etc.

    ! set module scope variables pointing to grid contents for ease of use
    call x%set_periodic_norefloc(flagperiodic)
    call set_gridflag(x%gridflag)
    call bind_grav_ptrs(x%g1,x%g2,x%g3)

    ! FIXME: do we need this???
    !> Make sure we have a sensible x2,3 decomposition of grid
    !> and that parameters aren't impossible
    !if(mpi_cfg%myid == 0) call grid_check(x)    
  end subroutine create_local_grid
end module grid
