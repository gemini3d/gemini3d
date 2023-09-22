module multifluid_mpi

use phys_consts, only: wp
use mpimod, only: mpi_cfg, tag=>gemini_mpi, halo

implicit none (type, external)
private
public :: halo_allparams, halo_fluidvars

contains
  !> halo all parameters that are solve time-dependently in preparation for advection
  subroutine halo_allparams(ns,rhovs1,rhoes,flagperiodic)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes
    logical, intent(in) :: flagperiodic
  
    call halo(ns,2,tag%ns,flagperiodic)
    call halo(rhovs1,2,tag%vs1,flagperiodic)
    call halo(rhoes,2,tag%Ts,flagperiodic)
  end subroutine halo_allparams


  !> halo all variables, including the perpendicular velocities
   subroutine halo_fluidvars(ns,rhovs1,rhoes,vs2,vs3,flagperiodic)
     real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes,vs2,vs3
     logical, intent(in) :: flagperiodic
     integer :: isp,lsp
     real(wp), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2) :: param

     call halo_allparams(ns,rhovs1,rhoes,flagperiodic)
     lsp=size(ns,4)
     do isp=1,lsp
       param=vs2(:,:,:,isp)
       call halo(param,2,tag%vs2BC,flagperiodic)
       vs2(:,:,:,isp)=param
     end do
     do isp=1,lsp
       param=vs3(:,:,:,isp)
       call halo(param,2,tag%vs3BC,flagperiodic)
       vs3(:,:,:,isp)=param
     end do
   end subroutine halo_fluidvars
end module multifluid_mpi
