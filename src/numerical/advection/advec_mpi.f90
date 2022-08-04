module advec_mpi

use phys_consts, only: wp
  !! do not import grid sizes in case we want do subgrid advection...
use mpimod, only: mpi_cfg, halo, tag=>gemini_mpi

implicit none (type, external)
private
public :: halo_interface_vels_allspec
contains
  !> Perform haloing needed to ghost-fill so cell interface vels (single species) can be computed across the grid
  subroutine halo_interface_vels(isp,isperiodic,vs2,vs3)
    integer, intent(in) :: isp
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
    real(wp), dimension(-1:size(vs3,1)-2,-1:size(vs3,2)-2,-1:size(vs3,3)-2) :: param

    !> NEED TO ALSO PASS THE X2 VELOCITIES SO WE CAN COMPUTE INTERFACE VALUES
    param=vs2(:,:,:,isp)
    call halo(param,1,tag%vs2BC,isperiodic)
    !! we only need one ghost cell to compute interface velocities
    vs2(:,:,:,isp)=param

    !> PASS X3 VELOCITY BOUNDARY CONDITIONS WITH GENERIC HALOING ROUTINES
    param=vs3(:,:,:,isp)
    call halo(param,1,tag%vs3BC,isperiodic)
    !! we only need one ghost cell to compute interface velocities
    vs3(:,:,:,isp)=param
  end subroutine halo_interface_vels


  !> Repeat haloing operations for all species.
  subroutine halo_interface_vels_allspec(isperiodic,vs2,vs3,lsp)
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
    integer, intent(in) :: lsp
    integer :: isp

    if (lsp>size(vs2,4)) error stop 'number of haloed species must be less than or equal to total species number'
    do isp=1,lsp
      call halo_interface_vels(isp,isperiodic,vs2,vs3)
    end do
  end subroutine halo_interface_vels_allspec
end module advec_mpi
