module neutral_perturbations

use phys_consts, only: wp, lnchem, pi, Re, debug
use config, only: gemini_cfg
use meshobj, only: curvmesh
use mpimod, only: mpi_cfg
use neutraldataobj, only: neutraldata
!!! FIXME: this should be conditioned on actually compiling with libgemini-mpi
use neutraldata3Dobj_mpi, only: neutraldata3D
!!!
use neutraldata2Daxisymmobj, only: neutraldata2Daxisymm
use neutraldata2Dcartobj, only: neutraldata2Dcart
use neutral, only: nnmsis,Tnmsis,vn1base,vn2base,vn3base

implicit none (type, external)

! flag to check whether to apply neutral perturbations
logical :: flagneuperturb=.false.

private
public :: init_neutralperturb,neutral_update,neutral_perturb,neutral_denstemp_update,neutral_wind_update,clear_dneu

contains
  !> initialize/allocate neutral perturbation data object
  subroutine init_neutralperturb(cfg,x,dt,ymd,UTsec,atmosperturb)
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
        allocate(neutraldata3D::atmosperturb)
      !!!
      case default
        error stop 'non-standard neutral interpolation type chosen in config.nml...'
      end select

      ! call object init procedure
      call atmosperturb%init(cfg,cfg%sourcedir,x,dt,cfg%dtneu,ymd,UTsec)
    end if
  end subroutine init_neutralperturb


 !> update neutral perturbations and add to main neutral arrays
  subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3,atmosperturb)
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dt,dtneu
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd
    !! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(inout) :: x
    !! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    !! intent(out)
    !! neutral params interpolated to plasma grid at requested time
    real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
    !! intent(out)
    class(neutraldata), intent(inout) :: atmosperturb

    ! advance object state
    call atmosperturb%update(cfg,dt,t,x,ymd,UTsec)

    !Add interpolated perturbations to module reference atmosphere arrays
    call neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)
  end subroutine neutral_perturb


  !> update density, temperature, and winds
  subroutine neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid,atmosperturb)
    !! adds stored base and perturbation neutral atmospheric parameters
    !!  these are module-scope parameters so not needed as input
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: Tn
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
    !! intent(out)
    real(wp) :: v2grid,v3grid
    class(neutraldata), intent(inout) :: atmosperturb

    call neutral_denstemp_update(nn,Tn,atmosperturb)
    call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid,atmosperturb)
  end subroutine neutral_update


  !> Adds stored base (viz. background) and perturbation neutral atmospheric density
  !   This does not use any of the existing data in arrays, but is declared inout to avoid potential
  !   issues with deallocation/reallocation.  This procedure should be used when you have both neutral
  !   background and perturbations and an update needs to be done. 
  subroutine neutral_denstemp_update(nn,Tn,atmosperturb)
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    real(wp), dimension(:,:,:), intent(inout) :: Tn
    class(neutraldata), intent(in) :: atmosperturb

    !> background neutral parameters
    nn=nnmsis
    Tn=Tnmsis

    !> add perturbations, if used
    if (flagneuperturb) then
      nn(:,:,:,1)=nn(:,:,:,1)+atmosperturb%dnOinow
      nn(:,:,:,2)=nn(:,:,:,2)+atmosperturb%dnN2inow
      nn(:,:,:,3)=nn(:,:,:,3)+atmosperturb%dnO2inow
      nn(:,:,:,1)=max(nn(:,:,:,1),1._wp)
      nn(:,:,:,2)=max(nn(:,:,:,2),1._wp)
      nn(:,:,:,3)=max(nn(:,:,:,3),1._wp)
      !! note we are not adjusting derived densities like NO since it's not clear how they may be related to major
      !! species perturbations.

      Tn=Tn+atmosperturb%dTninow
      Tn=max(Tn,50._wp)
    end if
  end subroutine neutral_denstemp_update


  !> update wind variables with background and perturbation quantities, note that this does not use any of the
  !    existing data in vn; but we still use intent(inout) to avoid weirdness with allocatable arrays.  This procedure
  !    should only be used when one has both a background and perturbation.  
  subroutine neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid,atmosperturb)
    real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
    real(wp), intent(in) :: v2grid,v3grid
    class(neutraldata), intent(in) :: atmosperturb

    !> background neutral parameters
    vn1=vn1base
    vn2=vn2base
    vn3=vn3base

    !> perturbations, if used
    if (flagneuperturb) then
      vn1=vn1+atmosperturb%dvn1inow
      vn2=vn2+atmosperturb%dvn2inow
      vn3=vn3+atmosperturb%dvn3inow
    end if

    !> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
    vn2=vn2-v2grid
    vn3=vn3-v3grid
  end subroutine neutral_wind_update


  !> deallocate neutral data object
  subroutine clear_dneu(atmosperturb)
    class(neutraldata), pointer, intent(inout) :: atmosperturb
 
    if (associated(atmosperturb)) deallocate(atmosperturb)
  end subroutine clear_dneu
end module neutral_perturbations
