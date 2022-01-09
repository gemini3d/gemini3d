!> non-mpi related grid quantities are stored here for use by other modules.  This could arguably be an 
!    object...
module grid

use phys_consts, only: wp

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
             set_total_grid_sizes,set_subgrid_sizes,set_gridflag,bind_grav_ptrs

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
end module grid
