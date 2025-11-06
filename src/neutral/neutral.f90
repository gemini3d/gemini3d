module neutral

!> This module performs and organizes top-level operations on neutral data used in the GEMINI model.  
!    Calls to MSIS, HWM etc. for background states and file- or p4est-based inputs for perturbations
!    are handled in separate modules

use phys_consts, only: wp, lnchem, pi, Re, debug
use grid, only: lx1, lx2, lx3
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use neutraldataobj, only: neutraldata

! also links MSIS from vendor/msis00/

implicit none (type, external)
private
public :: neutral_info, neutral_info_dealloc, neutral_info_alloc, &
  rotate_geo2native, store_geo2native_projections, rotate_native2geo, &
  neutral_aggregate, neutral_wind_aggregate


!> type encapsulating information needed by neutral module; FIXME: arguably this should be a class with type-bound procedures.
type neutral_info
  !! aggregate atmospheric state (BG + perturbations)
  real(wp), dimension(:,:,:,:), allocatable :: nn           ! neutral density array
  real(wp), dimension(:,:,:), allocatable :: Tn             ! neutral temperature
  real(wp), dimension(:,:,:), allocatable :: vn1,vn2,vn3    ! neutral velocities in model native components

  !! base msis state
  real(wp), dimension(:,:,:,:), allocatable :: nnBG
  real(wp), dimension(:,:,:), allocatable :: TnBG
  real(wp), dimension(:,:,:), allocatable :: vn1BG,vn2BG,vn3BG

  !! projection factors for converting vectors mag->geo; e.g. defining rotation matrix from geographic coords into
  logical :: flagprojections=.false.
  real(wp), dimension(:,:,:), allocatable :: proj_ealt_e1,proj_ealt_e2,proj_ealt_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eglat_e1,proj_eglat_e2,proj_eglat_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eglon_e1,proj_eglon_e2,proj_eglon_e3
end type neutral_info


contains
  !> update density, temperature, and winds in aggregate variable (BG + perturb)
  subroutine neutral_aggregate(v2grid,v3grid,atmos,atmosperturb,flagBGonly)
    !! adds stored base and perturbation neutral atmospheric parameters
    !!  these are module-scope parameters so not needed as input
    real(wp) :: v2grid,v3grid
    type(neutral_info), intent(inout) :: atmos
    class(neutraldata), pointer, intent(inout) :: atmosperturb
    logical, intent(in), optional :: flagBGonly

    if (present(flagBGonly)) then
      call neutral_denstemp_aggregate(atmos,atmosperturb,flagBGonly)
      call neutral_wind_aggregate(v2grid,v3grid,atmos,atmosperturb,flagBGonly)
    else
      call neutral_denstemp_aggregate(atmos,atmosperturb,.false.)
      call neutral_wind_aggregate(v2grid,v3grid,atmos,atmosperturb,.false.)
    end if
  end subroutine neutral_aggregate


  !> Adds stored base (viz. background) and perturbation neutral atmospheric density
  !   This does not use any of the existing data in arrays, but is declared inout to avoid potential
  !   issues with deallocation/reallocation.  This procedure should be used when you have both neutral
  !   background and perturbations and an update needs to be done.
  subroutine neutral_denstemp_aggregate(atmos,atmosperturb,flagBGonly)
    type(neutral_info), intent(inout) :: atmos
    class(neutraldata), pointer, intent(inout) :: atmosperturb
    logical, intent(in) :: flagBGonly

    !> background neutral parameters
    atmos%nn=atmos%nnBG
    atmos%Tn=atmos%TnBG

    !> add perturbations, if used
    if (associated(atmosperturb) .and. .not. flagBGonly) then
      atmos%nn(:,:,:,1)=atmos%nn(:,:,:,1)+atmosperturb%dnOinow
      atmos%nn(:,:,:,2)=atmos%nn(:,:,:,2)+atmosperturb%dnN2inow
      atmos%nn(:,:,:,3)=atmos%nn(:,:,:,3)+atmosperturb%dnO2inow
      atmos%nn(:,:,:,1)=max(atmos%nn(:,:,:,1),1._wp)
      atmos%nn(:,:,:,2)=max(atmos%nn(:,:,:,2),1._wp)
      atmos%nn(:,:,:,3)=max(atmos%nn(:,:,:,3),1._wp)
      !! note we are not adjusting derived densities like NO since it's not clear how they may be related to major
      !! species perturbations.

      atmos%Tn=atmos%Tn+atmosperturb%dTninow
      atmos%Tn=max(atmos%Tn,50._wp)
    end if
  end subroutine neutral_denstemp_aggregate


  !> update wind variables with background and perturbation quantities, note that this does not use any of the
  !    existing data in vn; but we still use intent(inout) to avoid weirdness with allocatable arrays.  This procedure
  !    should only be used when one has both a background and perturbation.
  subroutine neutral_wind_aggregate(v2grid,v3grid,atmos,atmosperturb,flagBGonly)
    real(wp), intent(in) :: v2grid,v3grid
    type(neutral_info), intent(inout) :: atmos
    class(neutraldata), pointer, intent(inout) :: atmosperturb
    logical, intent(in) :: flagBGonly

    !> background neutral parameters
    atmos%vn1=atmos%vn1BG
    atmos%vn2=atmos%vn2BG
    atmos%vn3=atmos%vn3BG

    !> perturbations, if used; already rotated into native coordinates
    if (associated(atmosperturb) .and. .not. flagBGonly) then
      atmos%vn1=atmos%vn1+atmosperturb%dvn1inow
      atmos%vn2=atmos%vn2+atmosperturb%dvn2inow
      atmos%vn3=atmos%vn3+atmosperturb%dvn3inow
    end if

    !> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
    !print*, '<><><><><><><><><> subtracting off:  ',v2grid,v3grid
    atmos%vn2=atmos%vn2-v2grid
    atmos%vn3=atmos%vn3-v3grid
  end subroutine neutral_wind_aggregate


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

    !> rotate vectors into model native coordinate system; check whether to store in base vs. perturb
    !    fields of the atmos derived type
    if (present(flagBG) .and. flagBG) then
      atmos%vn1BG=vnalt*atmos%proj_ealt_e1+vnglat*atmos%proj_eglat_e1+vnglon*atmos%proj_eglon_e1
      atmos%vn2BG=vnalt*atmos%proj_ealt_e2+vnglat*atmos%proj_eglat_e2+vnglon*atmos%proj_eglon_e2
      atmos%vn3BG=vnalt*atmos%proj_ealt_e3+vnglat*atmos%proj_eglat_e3+vnglon*atmos%proj_eglon_e3
    else
      atmos%vn1=vnalt*atmos%proj_ealt_e1+vnglat*atmos%proj_eglat_e1+vnglon*atmos%proj_eglon_e1
      atmos%vn2=vnalt*atmos%proj_ealt_e2+vnglat*atmos%proj_eglat_e2+vnglon*atmos%proj_eglon_e2
      atmos%vn3=vnalt*atmos%proj_ealt_e3+vnglat*atmos%proj_eglat_e3+vnglon*atmos%proj_eglon_e3
    end if
  end subroutine rotate_geo2native


  !> rotate winds from model native coordinate system (x1,x2,x3) to geographic coordinates.
  !    In this function we do not necessarily want to assume a location for input or output
  !    wind components; they should instead be provided as array inputs.  
  subroutine rotate_native2geo(vn1,vn2,vn3,vnalt,vnglat,vnglon,x,atmos)
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: vn1,vn2,vn3
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: vnalt,vnglat,vnglon
    class(curvmesh), intent(in) :: x
    type(neutral_info), intent(inout) :: atmos
    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3,3) :: ealt,eglat,eglon
    integer :: ix1,ix2,ix3,lx1,lx2,lx3

    lx1=x%lx1; lx2=x%lx2; lx3=x%lx3

    !> if first time called then allocate space for projections and compute
    if (.not. atmos%flagprojections) then
      call x%calc_unitvec_geo(ealt,eglon,eglat)
      call store_geo2native_projections(x,ealt,eglon,eglat,atmos)
      atmos%flagprojections=.true.
    end if

    !> rotate vectors into model native coordinate system; check whether to store in base vs. perturb
    !    fields of the atmos derived type.  Note tha the projection fields have size (lx1,lx2,lx3)
    !    whereas the size of the input and output arrays includes ghost cells...
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          vnalt(ix1,ix2,ix3)=vn1(ix1,ix2,ix3)*atmos%proj_ealt_e1(ix1,ix2,ix3) + &
                               vn2(ix1,ix2,ix3)*atmos%proj_ealt_e2(ix1,ix2,ix3) + &
                               vn3(ix1,ix2,ix3)*atmos%proj_ealt_e3(ix1,ix2,ix3)
          vnglon(ix1,ix2,ix3)=vn1(ix1,ix2,ix3)*atmos%proj_eglon_e1(ix1,ix2,ix3) + &
                               vn2(ix1,ix2,ix3)*atmos%proj_eglon_e2(ix1,ix2,ix3) + &
                               vn3(ix1,ix2,ix3)*atmos%proj_eglon_e3(ix1,ix2,ix3)
          vnglat(ix1,ix2,ix3)=vn1(ix1,ix2,ix3)*atmos%proj_eglat_e1(ix1,ix2,ix3) + &
                               vn2(ix1,ix2,ix3)*atmos%proj_eglat_e2(ix1,ix2,ix3) + &
                               vn3(ix1,ix2,ix3)*atmos%proj_eglat_e3(ix1,ix2,ix3)
        end do
      end do
    end do
  end subroutine rotate_native2geo


  !> compute projections for rotating winds geographic to native coordinate system
  subroutine store_geo2native_projections(x,ealt,eglon,eglat,atmos,rotmat)
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(in) :: ealt,eglon,eglat
    type(neutral_info), intent(inout) :: atmos
    real(wp), dimension(:,:,:,:,:), intent(inout), optional :: rotmat    ! for debugging purposes
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
    allocate(atmos%nnBG(lx1,lx2,lx3,lnchem),atmos%TnBG(lx1,lx2,lx3), &
               atmos%vn1BG(lx1,lx2,lx3),atmos%vn2BG(lx1,lx2,lx3),atmos%vn3BG(lx1,lx2,lx3))

    ! start everyone out at zero
    atmos%nn=0
    atmos%Tn=0
    atmos%vn1=0
    atmos%vn2=0
    atmos%vn3=0
    atmos%nnBG = 0
    atmos%TnBG = 0
    atmos%vn1BG = 0
    atmos%vn2BG = 0
    atmos%vn3BG = 0

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
    deallocate(atmos%nnBG,atmos%TnBG,atmos%vn1BG,atmos%vn2BG,atmos%vn3BG)

    ! rotations of neutral winds
    deallocate(atmos%proj_ealt_e1,atmos%proj_eglat_e1,atmos%proj_eglon_e1)
    deallocate(atmos%proj_ealt_e2,atmos%proj_eglat_e2,atmos%proj_eglon_e2)
    deallocate(atmos%proj_ealt_e3,atmos%proj_eglat_e3,atmos%proj_eglon_e3)
    atmos%flagprojections=.false.
  end subroutine neutral_info_dealloc
end module neutral
