module neutral

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: x, lx1, lx2, lx3
use meshobj, only: curvmesh
use timeutils, only : find_lastdate
use config, only: gemini_cfg,cfg

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: nnmsis,Tnmsis, neutral_atmos, make_neuBG, clear_neuBG, init_neutralBG, &
  neutral_winds, rotate_geo2native, store_geo2native_projections, neutralBG_denstemp, neutralBG_wind, &
  vn1base,vn2base,vn3base

interface !< atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,msis_version)
    integer, intent(in) :: ymd(3), msis_version
    real(wp), intent(in) :: UTsecd
    real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
    real(wp), intent(in) :: activ(3)
  end subroutine neutral_atmos
end interface
interface !< wind.f90
  module subroutine neutral_winds(ymd, UTsec, Ap, x)
    integer, intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec, Ap
    class(curvmesh), intent(in) :: x
  end subroutine neutral_winds
end interface

!! BASE MSIS ATMOSPHERIC STATE ON WHICH TO APPLY PERTURBATIONS
real(wp), dimension(:,:,:,:), allocatable, protected :: nnmsis
real(wp), dimension(:,:,:), allocatable, protected :: Tnmsis
real(wp), dimension(:,:,:), allocatable, protected :: vn1base,vn2base,vn3base

!! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3

contains
  !> initializes neutral atmosphere by:
  !    1)  allocating storage space
  !    2)  establishing initial background for density, temperature, and winds
  !! Arguably this should be called for init and for neutral updates, except for the allocation part...
  subroutine init_neutralBG(dt,t,ymd,UTsec,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
    real(wp), intent(in) :: dt,t
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
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
    call make_neuBG()

    !! call msis to get an initial neutral background atmosphere
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call cpu_time(tstart)
    call neutral_atmos(cfg%ymd0,cfg%UTsec0,x%glat,x%glon,x%alt,cfg%activ,cfg%msis_version)
    call neutralBG_denstemp(nn,Tn)
    !if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if

    !> Horizontal wind model initialization/background
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call cpu_time(tstart)
    call neutral_winds(cfg%ymd0, cfg%UTsec0, Ap=cfg%activ(3), x=x)
    call neutralBG_wind(vn1,vn2,vn3,v2grid,v3grid)
    !! we sum the horizontal wind with the background state vector
    !! if HWM14 is disabled, neutral_winds returns the background state vector unmodified
    !if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral winds (from HWM) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if
  end subroutine init_neutralBG


  !>  Sets the neutral density and temperature variables to background values (usually from MSIS).  Do not use this procedure
  !     when using neutral perturbations; there is a sister procedure in that module for doing updates in that case.
  subroutine neutralBG_denstemp(nn,Tn)
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    real(wp), dimension(:,:,:), intent(inout) :: Tn

    !> background neutral parameters
    nn=nnmsis
    Tn=Tnmsis
  end subroutine neutralBG_denstemp


  !>  Sets neutral winds to background values (usually from MSIS).  If running with neutral perturbations use the sister
  !     procedure in neutral_perturb instead of this one.  
  subroutine neutralBG_wind(vn1,vn2,vn3,v2grid,v3grid)
    real(wp), dimension(:,:,:), intent(inout) :: vn1,vn2,vn3
    real(wp), intent(in) :: v2grid,v3grid

    !> background neutral parameters
    vn1=vn1base
    vn2=vn2base
    vn3=vn3base

    !> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
    vn2=vn2-v2grid
    vn3=vn3-v3grid
  end subroutine neutralBG_wind
 

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


  subroutine make_neuBG()
    !allocate and compute plasma grid z,rho locations and space to save neutral perturbation variables and projection factors
    allocate(nnmsis(lx1,lx2,lx3,lnchem),Tnmsis(lx1,lx2,lx3),vn1base(lx1,lx2,lx3),vn2base(lx1,lx2,lx3),vn3base(lx1,lx2,lx3))

    !start everyone out at zero
    nnmsis = 0
    Tnmsis = 0
    vn1base = 0
    vn2base = 0
    vn3base = 0
  end subroutine make_neuBG


  subroutine clear_neuBG
    ! stuff allocated at beginning of program
    deallocate(nnmsis,Tnmsis,vn1base,vn2base,vn3base)

    ! rotations of neutral winds
    if (allocated(proj_ealt_e1)) then
      deallocate(proj_ealt_e1,proj_eglat_e1,proj_eglon_e1)
      deallocate(proj_ealt_e2,proj_eglat_e2,proj_eglon_e2)
      deallocate(proj_ealt_e3,proj_eglat_e3,proj_eglon_e3)
    end if
  end subroutine clear_neuBG
end module neutral
