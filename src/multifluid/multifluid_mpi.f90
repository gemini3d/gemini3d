module multifluid_mpi

use phys_consts, only: wp
use mpimod, only: mpi_cfg, tag=>gemini_mpi, halo

implicit none (type, external)
private
public :: halo_allparams

contains

!> halo all parameters in preparation for advection
subroutine halo_allparams(ns,rhovs1,rhoes,flagperiodic)
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,rhoes
  logical, intent(in) :: flagperiodic

  call halo(ns,2,tag%ns,flagperiodic)
  call halo(rhovs1,2,tag%vs1,flagperiodic)
  call halo(rhoes,2,tag%Ts,flagperiodic)
end subroutine halo_allparams

end module multifluid_mpi
