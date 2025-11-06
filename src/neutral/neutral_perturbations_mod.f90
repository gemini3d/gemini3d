module neutral_perturbations

use phys_consts, only: wp, lnchem, pi, Re, debug
use gemini3d_config, only: gemini_cfg
use meshobj, only: curvmesh
use mpimod, only: mpi_cfg
use neutraldataobj, only: neutraldata
!!! FIXME: this should be conditioned on actually compiling with libgemini-mpi
!use neutraldata3Dobj_mpi, only: neutraldata3D
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use neutraldata3Dobj_fclaw_axisymm, only: neutraldata3D_fclaw_axisymm
use neutraldata3Dobj_fclaw_3Dx, only: neutraldata3D_fclaw_3Dx
use neutraldata3Dobj_geog_mpi, only: neutraldata3D_geog
use neutraldata3Dobj_geom_mpi, only: neutraldata3D_geom
!!!
use neutraldata2Daxisymmobj, only: neutraldata2Daxisymm
use neutraldata2Dcartobj, only: neutraldata2Dcart
use neutral, only: neutral_info

implicit none (type, external)

! flag to check whether to apply neutral perturbations
logical :: flagneuperturb=.false.

private
public :: init_neutral_perturb,neutral_perturb,clear_neutral_perturb

contains
  !> initialize/allocate neutral perturbation data object
  subroutine init_neutral_perturb(cfg,x,dt,ymd,UTsec,atmosperturb)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dt
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(neutraldata), pointer, intent(inout) :: atmosperturb

    !! perform an initialization for the perturbation quantities
    if (cfg%flagdneu==1) then
      ! set flag denoted neutral perturbations
      flagneuperturb=.true.

      ! allocate correct type, FIXME: eventuallly no shunt to 3D
      select case (cfg%interptype)
      case (0)
        allocate(neutraldata2Dcart::atmosperturb)
      case (1)
        allocate(neutraldata2Daxisymm::atmosperturb)
      !!! FIXME: conditioned on compiling with libgemini-mpi
      case (3)
        allocate(neutraldata3D_geom::atmosperturb)
      case (4)
        allocate(neutraldata3D_geog::atmosperturb)
      case (5)     ! we assume forestGEMINI will control things, this will be axisymmetric MAGIC-forest data
        allocate(neutraldata3D_fclaw_axisymm::atmosperturb)
      case (6)     ! we assume forestGEMINI will control things, this will be 3Dx MAGIC-forest data
        allocate(neutraldata3D_fclaw_3Dx::atmosperturb)
      case default
        print*, 'non-standard neutral interpolation type chosen in config.nml:  ',cfg%interptype
        error stop
      end select

      ! call object init procedure
      call atmosperturb%init(cfg,cfg%sourcedir,x,dt,cfg%dtneu,ymd,UTsec)
    end if
  end subroutine init_neutral_perturb


 !> update neutral perturbations and add to main neutral arrays.  It is assumed that the main GEMINI app will control
 !    when and where this gets called.
  subroutine neutral_perturb(cfg,dt,t,ymd,UTsec,x,v2grid,v3grid,atmos,atmosperturb)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dt
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    !! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(inout) :: x
    !! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
    real(wp), intent(in) :: v2grid,v3grid
    type(neutral_info), intent(inout) :: atmos
    class(neutraldata), pointer, intent(inout) :: atmosperturb

    ! advance object state
    call atmosperturb%update(cfg,dt,t,x,ymd,UTsec)

    !Add interpolated perturbations to module reference atmosphere arrays -- move to libgemini_mpi
    !call neutral_aggregate(v2grid,v3grid,atmos,atmosperturb)
  end subroutine neutral_perturb


  !> deallocate neutral data object
  subroutine clear_neutral_perturb(atmosperturb)
    class(neutraldata), pointer, intent(inout) :: atmosperturb

    if (associated(atmosperturb)) deallocate(atmosperturb)
  end subroutine clear_neutral_perturb
end module neutral_perturbations
