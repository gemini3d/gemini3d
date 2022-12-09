module neutral

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: lx1, lx2, lx3
use meshobj, only: curvmesh
use timeutils, only : find_lastdate
use gemini3d_config, only: gemini_cfg

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: neutral_atmos, init_neutralBG, neutral_info, neutral_info_dealloc, neutral_info_alloc, &
  neutral_winds, rotate_geo2native, store_geo2native_projections, neutralBG_denstemp, neutralBG_wind


!> type encapsulating information needed by neutral module; FIXME: arguably this should be a class with type-bound procedures.
type neutral_info
  !! total//overall atmospheric state (BG + perturbations)
  real(wp), dimension(:,:,:,:), allocatable :: nn    !! neutral density array
  real(wp), dimension(:,:,:), allocatable :: Tn
  real(wp), dimension(:,:,:), allocatable :: vn1,vn2,vn3    !! neutral temperature and velocities

  !! base msis state
  real(wp), dimension(:,:,:,:), allocatable :: nnmsis
  real(wp), dimension(:,:,:), allocatable :: Tnmsis
  real(wp), dimension(:,:,:), allocatable :: vn1base,vn2base,vn3base

  !! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
  logical :: flagprojections=.false.
  real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3
end type neutral_info


interface !< atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,msis_version,atmos)
    integer, intent(in) :: ymd(3), msis_version
    real(wp), intent(in) :: UTsecd
    real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
    real(wp), intent(in) :: activ(3)
    type(neutral_info), intent(inout) :: atmos
  end subroutine neutral_atmos
end interface
interface !< wind.f90
  module subroutine neutral_winds(ymd, UTsec, Ap, x, atmos)
    integer, intent(in) :: ymd(3)
    real(wp), intent(in) :: UTsec, Ap
    class(curvmesh), intent(in) :: x
    type(neutral_info), intent(inout) :: atmos
  end subroutine neutral_winds
end interface

contains
  !> initializes neutral atmosphere by:
  !    1)  establishing initial background for density, temperature, and winds
  !! Arguably this should be called for init and for neutral updates, except for the allocation part...
  subroutine init_neutralBG(cfg,ymd,UTsec,x,v2grid,v3grid,atmos)

    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(inout) :: x    ! unit vecs may be deallocated after first setup
    real(wp), intent(in) :: v2grid,v3grid
    type(neutral_info), intent(inout) :: atmos
    real(wp) :: tstart,tfin

    !! allocation neutral module scope variables so there is space to store all the file input and do interpolations
    !call make_neuBG()

    !! call msis to get an initial neutral background atmosphere
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call cpu_time(tstart)
    call neutral_atmos(cfg%ymd0,cfg%UTsec0,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3), &
                         cfg%activ,cfg%msis_version,atmos)
    call neutralBG_denstemp(atmos)
    !if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if

    !> Horizontal wind model initialization/background
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    call cpu_time(tstart)
    call neutral_winds(cfg%ymd0, cfg%UTsec0, Ap=cfg%activ(3), x=x, atmos=atmos)
    call neutralBG_wind(atmos,v2grid,v3grid)
    !! we sum the horizontal wind with the background state vector
    !! if HWM14 is disabled, neutral_winds returns the background state vector unmodified
    !if (mpi_cfg%myid == 0) then
      call cpu_time(tfin)
      print *, 'Initial neutral winds (from HWM) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if
  end subroutine init_neutralBG


  !>  Sets the neutral density and temperature variables to background values (usually from MSIS).  Do not use this procedure
  !     when using neutral perturbations; there is a sister procedure in that module for doing updates in that case.
  subroutine neutralBG_denstemp(atmos)
    type(neutral_info), intent(inout) :: atmos

    !> background neutral parameters
    atmos%nn=atmos%nnmsis
    atmos%Tn=atmos%Tnmsis
  end subroutine neutralBG_denstemp


  !>  Sets neutral winds to background values (usually from MSIS).  If running with neutral perturbations use the sister
  !     procedure in neutral_perturb instead of this one.
  subroutine neutralBG_wind(atmos,v2grid,v3grid)
    type(neutral_info), intent(inout) :: atmos
    real(wp), intent(in) :: v2grid,v3grid

    !> background neutral parameters
    atmos%vn1=atmos%vn1base
    atmos%vn2=atmos%vn2base
    atmos%vn3=atmos%vn3base

    !> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
    atmos%vn2=atmos%vn2-v2grid
    atmos%vn3=atmos%vn3-v3grid
  end subroutine neutralBG_wind


  !> rotate winds from geographic to model native coordinate system (x1,x2,x3)
  subroutine rotate_geo2native(vnalt,vnglat,vnglon,x,atmos,flagBG)
    real(wp), dimension(:,:,:), intent(in) :: vnalt,vnglat,vnglon
    class(curvmesh), intent(in) :: x
    type(neutral_info), intent(inout) :: atmos
    logical, intent(in), optional :: flagBG
    real(wp), dimension(1:size(vnalt,1),1:size(vnalt,2),1:size(vnalt,3),3) :: ealt,eglat,eglon

    !> if first time called then allocate space for projections and compute
    if (.not. atmos%flagprojections) then
      call x%calc_unitvec_geo(ealt,eglon,eglat)
      call store_geo2native_projections(x,ealt,eglon,eglat,atmos)
      atmos%flagprojections=.true.
    end if

    !> rotate vectors into model native coordinate system
    if (present(flagBG) .and. flagBG) then
      atmos%vn1base=vnalt*atmos%proj_ealt_e1+vnglat*atmos%proj_eglat_e1+vnglon*atmos%proj_eglon_e1
      atmos%vn2base=vnalt*atmos%proj_ealt_e2+vnglat*atmos%proj_eglat_e2+vnglon*atmos%proj_eglon_e2
      atmos%vn3base=vnalt*atmos%proj_ealt_e3+vnglat*atmos%proj_eglat_e3+vnglon*atmos%proj_eglon_e3
    else
      atmos%vn1=vnalt*atmos%proj_ealt_e1+vnglat*atmos%proj_eglat_e1+vnglon*atmos%proj_eglon_e1
      atmos%vn2=vnalt*atmos%proj_ealt_e2+vnglat*atmos%proj_eglat_e2+vnglon*atmos%proj_eglon_e2
      atmos%vn3=vnalt*atmos%proj_ealt_e3+vnglat*atmos%proj_eglat_e3+vnglon*atmos%proj_eglon_e3
    end if
  end subroutine rotate_geo2native


  !> compute projections for rotating winds geographic to native coordinate system
  subroutine store_geo2native_projections(x,ealt,eglon,eglat,atmos,rotmat)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(in) :: ealt,eglon,eglat
    type(neutral_info), intent(inout) :: atmos
    real(wp), dimension(:,:,:,:,:), intent(out), optional :: rotmat    ! for debugging purposes
    integer :: ix1,ix2,ix3,lx1,lx2,lx3

    !! allocate module-scope space for the projection factors
    lx1=size(ealt,1); lx2=size(ealt,2); lx3=size(ealt,3);

    !! compute projections (dot products of unit vectors)
    atmos%proj_ealt_e1=sum(ealt*x%e1,4)
    atmos%proj_eglat_e1=sum(eglat*x%e1,4)
    atmos%proj_eglon_e1=sum(eglon*x%e1,4)
    atmos%proj_ealt_e2=sum(ealt*x%e2,4)
    atmos%proj_eglat_e2=sum(eglat*x%e2,4)
    atmos%proj_eglon_e2=sum(eglon*x%e2,4)
    atmos%proj_ealt_e3=sum(ealt*x%e3,4)
    atmos%proj_eglat_e3=sum(eglat*x%e3,4)
    atmos%proj_eglon_e3=sum(eglon*x%e3,4)

    !! store the rotation matrix to convert geo to native if the user wants it
    if (present(rotmat)) then
      do ix3=1,lx3
        do ix2=1,lx2
          do ix1=1,lx1
            rotmat(1,1:3,ix1,ix2,ix3)=[atmos%proj_ealt_e1(ix1,ix2,ix3),atmos%proj_eglat_e1(ix1,ix2,ix3), &
                                         atmos%proj_eglon_e1(ix1,ix2,ix3)]
            rotmat(2,1:3,ix1,ix2,ix3)=[atmos%proj_ealt_e2(ix1,ix2,ix3),atmos%proj_eglat_e2(ix1,ix2,ix3), &
                                         atmos%proj_eglon_e2(ix1,ix2,ix3)]
            rotmat(3,1:3,ix1,ix2,ix3)=[atmos%proj_ealt_e3(ix1,ix2,ix3),atmos%proj_eglat_e3(ix1,ix2,ix3), &
                                         atmos%proj_eglon_e3(ix1,ix2,ix3)]
          end do
        end do
      end do
    end if
  end subroutine store_geo2native_projections


  !> allocate space for data used in neutral module
  subroutine neutral_info_alloc(atmos)
    type(neutral_info), intent(inout) :: atmos

    ! allocate space for overall neutral state
    allocate(atmos%nn(lx1,lx2,lx3,lnchem),atmos%Tn(lx1,lx2,lx3),atmos%vn1(lx1,lx2,lx3), &
               atmos%vn2(lx1,lx2,lx3),atmos%vn3(lx1,lx2,lx3))

    ! allocate space for background atmospheric state
    allocate(atmos%nnmsis(lx1,lx2,lx3,lnchem),atmos%Tnmsis(lx1,lx2,lx3), &
               atmos%vn1base(lx1,lx2,lx3),atmos%vn2base(lx1,lx2,lx3),atmos%vn3base(lx1,lx2,lx3))

    ! start everyone out at zero
    atmos%nn=0
    atmos%Tn=0
    atmos%vn1=0
    atmos%vn2=0
    atmos%vn3=0
    atmos%nnmsis = 0
    atmos%Tnmsis = 0
    atmos%vn1base = 0
    atmos%vn2base = 0
    atmos%vn3base = 0

    ! projection factors for geographic rotations
    allocate(atmos%proj_ealt_e1(lx1,lx2,lx3),atmos%proj_eglat_e1(lx1,lx2,lx3),atmos%proj_eglon_e1(lx1,lx2,lx3))
    allocate(atmos%proj_ealt_e2(lx1,lx2,lx3),atmos%proj_eglat_e2(lx1,lx2,lx3),atmos%proj_eglon_e2(lx1,lx2,lx3))
    allocate(atmos%proj_ealt_e3(lx1,lx2,lx3),atmos%proj_eglat_e3(lx1,lx2,lx3),atmos%proj_eglon_e3(lx1,lx2,lx3))
  end subroutine neutral_info_alloc


  !> deallocate arrays in neutral_info type
  subroutine neutral_info_dealloc(atmos)
    type(neutral_info), intent(inout) :: atmos

    ! overall neutral data
    deallocate(atmos%nn,atmos%Tn,atmos%vn1,atmos%vn2,atmos%vn3)

    ! msis data
    deallocate(atmos%nnmsis,atmos%Tnmsis,atmos%vn1base,atmos%vn2base,atmos%vn3base)

    ! rotations of neutral winds
    deallocate(atmos%proj_ealt_e1,atmos%proj_eglat_e1,atmos%proj_eglon_e1)
    deallocate(atmos%proj_ealt_e2,atmos%proj_eglat_e2,atmos%proj_eglon_e2)
    deallocate(atmos%proj_ealt_e3,atmos%proj_eglat_e3,atmos%proj_eglon_e3)
    atmos%flagprojections=.false.
  end subroutine neutral_info_dealloc
end module neutral
