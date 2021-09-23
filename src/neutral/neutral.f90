module neutral

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: lx1, lx2, lx3
use meshobj, only: curvmesh
use timeutils, only : find_lastdate
use mpimod, only: mpi_cfg
use config, only: gemini_cfg
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata2Daxisymmobj, only: neutraldata2Daxisymm
use neutraldata2Dcartobj, only: neutraldata2Dcart

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: Tnmsis, neutral_atmos, make_dneu, clear_dneu, neutral_perturb, neutral_update, init_neutrals, &
  neutral_winds, rotate_geo2native, neutral_denstemp_update, neutral_wind_update, store_geo2native_projections


interface !< atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,nn,Tn,msis_version)
    integer, intent(in) :: ymd(3), msis_version
    real(wp), intent(in) :: UTsecd
    real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
    real(wp), intent(in) :: activ(3)
    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3),lnchem), intent(inout) :: nn
    !! intent(out)
    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(inout) :: Tn
    !! intent(out)
  end subroutine neutral_atmos
end interface

interface !< wind.f90
  module subroutine neutral_winds(ymd, UTsec, Ap, x, v2grid, v3grid, vn1, vn2, vn3)
    integer, intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec, Ap
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(1:size(x%alt,1),1:size(x%alt,2),1:size(x%alt,3)), intent(inout) :: vn1,vn2,vn3
  end subroutine neutral_winds
end interface

! flag to check whether to apply neutral perturbations
logical :: flagneuperturb=.false.

!! BASE MSIS ATMOSPHERIC STATE ON WHICH TO APPLY PERTURBATIONS
real(wp), dimension(:,:,:,:), allocatable, protected :: nnmsis
real(wp), dimension(:,:,:), allocatable, protected :: Tnmsis
real(wp), dimension(:,:,:), allocatable, protected :: vn1base,vn2base,vn3base

!! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3

!! new module variables for OO refactor, polymorphic perturbation object
class(neutraldata), allocatable :: atmosperturb

contains
  !> initializes neutral atmosphere by:
  !    1)  allocating storage space
  !    2)  establishing initial background for density, temperature, and winds
  !    3)  priming file input so that we have an initial perturbed state to start from (necessary for restart)
  subroutine init_neutrals(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
    real(wp), intent(in) :: dt,t
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(inout) :: x    ! unit vecs may be deallocated after first setup
    real(wp), intent(in) :: v2grid,v3grid
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: Tn
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
    !! intent(out)
  
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp
    real(wp) :: tstart,tfin
    
    !! allocation neutral module scope variables so there is space to store all the file input and do interpolations
    call make_dneu()
    
    !! call msis to get an initial neutral background atmosphere
    if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call neutral_atmos(cfg%ymd0,cfg%UTsec0,x%glat,x%glon,x%alt,cfg%activ,nn,Tn,cfg%msis_version)
    if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    end if
    
    !> Horizontal wind model initialization/background
    if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call neutral_winds(cfg%ymd0, cfg%UTsec0, Ap=cfg%activ(3), x=x, v2grid=v2grid,v3grid=v3grid,vn1=vn1, vn2=vn2, vn3=vn3)
    !! we sum the horizontal wind with the background state vector
    !! if HWM14 is disabled, neutral_winds returns the background state vector unmodified
    if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral winds (from HWM) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    end if
  
  
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
      case (3)
        allocate(neutraldata3D::atmosperturb)
      case default
        error stop 'non-standard neutral interpolation type chosen in config.nml...'
      end select
  
      ! call object init procedure
      call atmosperturb%init(cfg,cfg%sourcedir,x,dt,cfg%dtneu,ymd,UTsec)
    end if
  end subroutine init_neutrals


  !> update neutral perturbations and add to main neutral arrays
  module subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
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

    ! advance object state
    call atmosperturb%update(cfg,dt,t,x,ymd,UTsec)
  
    !Add interpolated perturbations to module reference atmosphere arrays
    call neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)
  end subroutine neutral_perturb
  
  
  !> update density, temperature, and winds
  subroutine neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)
    !! adds stored base and perturbation neutral atmospheric parameters
    !!  these are module-scope parameters so not needed as input
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    !! intent(out)  
    real(wp), dimension(:,:,:), intent(inout) :: Tn
    !! intent(out)
    real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
    !! intent(out)
    real(wp) :: v2grid,v3grid
  
    call neutral_denstemp_update(nn,Tn)
    call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
  end subroutine neutral_update
  
  
  !> Adds stored base (viz. background) and perturbation neutral atmospheric density
  subroutine neutral_denstemp_update(nn,Tn)
    real(wp), dimension(:,:,:,:), intent(out) :: nn
    real(wp), dimension(:,:,:), intent(out) :: Tn
    
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
  
  
  !> update wind variables with background and perturbation quantities
  subroutine neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
    real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
    real(wp) :: v2grid,v3grid
    
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
  
  
  !> rotate winds from geographic to model native coordinate system (x1,x2,x3)
  subroutine rotate_geo2native(vnalt,vnglat,vnglon,x,vn1,vn2,vn3)
    real(wp), dimension(:,:,:), intent(in) :: vnalt,vnglat,vnglon
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
    real(wp), dimension(1:size(vnalt,1),1:size(vnalt,2),1:size(vnalt,3),3) :: ealt,eglat,eglon
    integer :: lx1,lx2,lx3
  
    !> if first time called then allocate space for projections and compute
    if (.not. allocated(proj_ealt_e1)) then
      call x%calc_unitvec_geo(ealt,eglon,eglat)
      call store_geo2native_projections(x,ealt,eglon,eglat)
    end if
  
    !> rotate vectors into model native coordinate system
    vn1=vnalt*proj_ealt_e1+vnglat*proj_eglat_e1+vnglon*proj_eglon_e1
    vn2=vnalt*proj_ealt_e2+vnglat*proj_eglat_e2+vnglon*proj_eglon_e2
    vn3=vnalt*proj_ealt_e3+vnglat*proj_eglat_e3+vnglon*proj_eglon_e3
  end subroutine rotate_geo2native
  
  
  !> compute projections for rotating winds geographic to native coordinate system
  subroutine store_geo2native_projections(x,ealt,eglon,eglat,rotmat)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(in) :: ealt,eglon,eglat
    real(wp), dimension(:,:,:,:,:), intent(out), optional :: rotmat    ! for debugging purposes
    integer :: ix1,ix2,ix3,lx1,lx2,lx3
  
    !! allocate module-scope space for the projection factors
    lx1=size(ealt,1); lx2=size(ealt,2); lx3=size(ealt,3);
    allocate(proj_ealt_e1(lx1,lx2,lx3),proj_eglat_e1(lx1,lx2,lx3),proj_eglon_e1(lx1,lx2,lx3))
    allocate(proj_ealt_e2(lx1,lx2,lx3),proj_eglat_e2(lx1,lx2,lx3),proj_eglon_e2(lx1,lx2,lx3))
    allocate(proj_ealt_e3(lx1,lx2,lx3),proj_eglat_e3(lx1,lx2,lx3),proj_eglon_e3(lx1,lx2,lx3))
  
    !! compute projections (dot products of unit vectors)
    proj_ealt_e1=sum(ealt*x%e1,4)
    proj_eglat_e1=sum(eglat*x%e1,4)
    proj_eglon_e1=sum(eglon*x%e1,4)
    proj_ealt_e2=sum(ealt*x%e2,4)
    proj_eglat_e2=sum(eglat*x%e2,4)
    proj_eglon_e2=sum(eglon*x%e2,4)
    proj_ealt_e3=sum(ealt*x%e3,4)
    proj_eglat_e3=sum(eglat*x%e3,4)
    proj_eglon_e3=sum(eglon*x%e3,4)
  
    !! store the rotation matrix to convert geo to native if the user wants it
    if (present(rotmat)) then
      do ix3=1,lx3
        do ix2=1,lx2
          do ix1=1,lx1
            rotmat(1,1:3,ix1,ix2,ix3)=[proj_ealt_e1(ix1,ix2,ix3),proj_eglat_e1(ix1,ix2,ix3),proj_eglon_e1(ix1,ix2,ix3)]
            rotmat(2,1:3,ix1,ix2,ix3)=[proj_ealt_e2(ix1,ix2,ix3),proj_eglat_e2(ix1,ix2,ix3),proj_eglon_e2(ix1,ix2,ix3)]
            rotmat(3,1:3,ix1,ix2,ix3)=[proj_ealt_e3(ix1,ix2,ix3),proj_eglat_e3(ix1,ix2,ix3),proj_eglon_e3(ix1,ix2,ix3)]
          end do
        end do
      end do
    end if
  end subroutine store_geo2native_projections
  
  
  subroutine make_dneu()
    !allocate and compute plasma grid z,rho locations and space to save neutral perturbation variables and projection factors
    allocate(nnmsis(lx1,lx2,lx3,lnchem),Tnmsis(lx1,lx2,lx3),vn1base(lx1,lx2,lx3),vn2base(lx1,lx2,lx3),vn3base(lx1,lx2,lx3))
    
    !start everyone out at zero
    nnmsis = 0
    Tnmsis = 0
    vn1base = 0
    vn2base = 0
    vn3base = 0
  end subroutine make_dneu
  
  
  subroutine clear_dneu
    ! stuff allocated at beginning of program
    deallocate(nnmsis,Tnmsis,vn1base,vn2base,vn3base)
    
    ! rotations of neutral winds
    if (allocated(proj_ealt_e1)) then
      deallocate(proj_ealt_e1,proj_eglat_e1,proj_eglon_e1)
      deallocate(proj_ealt_e2,proj_eglat_e2,proj_eglon_e2)
      deallocate(proj_ealt_e3,proj_eglat_e3,proj_eglon_e3)
    end if
  end subroutine clear_dneu
end module neutral
