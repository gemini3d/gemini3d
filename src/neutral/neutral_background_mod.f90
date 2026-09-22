module neutral_background

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp
use neutral, only: rotate_geo2native, neutral_info
use neutraldataBGobj, only: neutraldataBG
use msis_interface, only : msisinit
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use grid, only: lx1,lx2,lx3

implicit none (type, external)
private
public :: init_neutral_background, neutral_background_fileinput, neutral_background_empirical, &
        clear_neutral_background_fileinput, msisinit_in

interface !< atmos.f90
  module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,msis_version,atmos)
    integer, intent(in) :: ymd(3), msis_version
    real(wp), intent(in) :: UTsecd
    real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
    real(wp), intent(in) :: activ(3)
    type(neutral_info), intent(inout) :: atmos
  end subroutine neutral_atmos
  module subroutine NO_calc(nn,Tn)
    real(wp), dimension(:,:,:,:), intent(inout) :: nn
    real(wp), dimension(:,:,:), intent(in) :: Tn
  end subroutine NO_calc
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
  !> Determine whether we are using MSIS/HWM or file-based background atmospheric information
  subroutine init_neutral_background(dt,cfg,ymd,UTsec,x,v2grid,v3grid,atmos,atmosbackground)
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(inout) :: x
    real(wp), intent(in) :: v2grid,v3grid
    type(neutral_info), intent(inout) :: atmos
    type(neutraldataBG), intent(inout) :: atmosbackground

    if (cfg%flagneutralBGfile==1) then
      print*, 'NOTE:  Neutral background to be taken from file input...'
      !allocate(neutralBGdata::atmos%atmosBGdata)    ! not a class pointer so is allocated a priori
      call atmosbackground%init(cfg,cfg%neutralBGdir,x,dt,cfg%dtneuBGfile,ymd,UTsec)
      call neutral_background_fileinput_copyout(x,atmos,atmosbackground)          ! need to move data from object to shared space used by others
    else
      print*, 'NOTE:  Neutral background to be taken from empirical models...'
      call msisinit_in(cfg)
      call neutral_background_empirical(cfg,ymd,UTsec,x,v2grid,v3grid,atmos)    ! this will actually store the data in shared space
    end if
  end subroutine init_neutral_background


  !> Neutral background values based on file input, then place into neutral_info struct
  subroutine neutral_background_fileinput(dtmodel,t,cfg,ymd,UTsec,x,atmos,atmosbackground)
    real(wp), intent(in) :: dtmodel
    real(wp), intent(in) :: t
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    type(neutral_info), intent(inout) :: atmos
    type(neutraldataBG), intent(inout) :: atmosbackground

    call atmosbackground%update(cfg,dtmodel,t,x,ymd,UTsec)

    call neutral_background_fileinput_copyout(x,atmos,atmosbackground)
  end subroutine neutral_background_fileinput


  subroutine neutral_background_fileinput_copyout(x,atmos,atmosbackground)
    class(curvmesh), intent(in) :: x
    type(neutral_info), intent(inout) :: atmos
    type(neutraldataBG), intent(inout) :: atmosbackground

    integer :: ix1,ix2,ix3,ineu,ix1ref1,ix1ref2,ix1apex
    integer :: lneu
    real(wp) :: Tnref,vn1ref,vn2ref,vn3ref,nrat

    ! place the background atmospheric density and temperature (drifts still need to be rotated)
    !   from object into arrays used by other modules
    atmos%nnBG(:,:,:,1:5)=atmosbackground%natminow(:,:,:,1:5)
    atmos%TnBG=atmosbackground%natminow(:,:,:,9)
    if (.not. all(ieee_is_finite(atmos%TnBG)) .or. any(atmos%TnBG<=0)) &
      error stop 'neutral background: temperatures must be finite and positive'
    call NO_calc(atmos%nnBG,atmos%TnBG)

    ! rotate winds to model native and place into wind background arrays in the atmos object
    call rotate_geo2native(atmosbackground%natminow(:,:,:,6), &
            atmosbackground%natminow(:,:,:,7), &
            atmosbackground%natminow(:,:,:,8), &
            x,atmos,.true.)

    lneu=size(atmos%nnBG,4)

    ! We need to extrapolate in altitude from the nearest hemisphere ionosphere in addition to dealing with
    !   issues that arise from below-ground cells.
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        ! for each field line we need to find the apex altitude and index (it possibly could be a different index
        !   on each field line???
        if (.not. all(ieee_is_finite(x%alt(1:lx1,ix2,ix3))) .or. &
            .not. ieee_is_finite(atmosbackground%altpmax)) error stop 'neutral background: nonfinite altitude'
        ix1apex=maxloc(x%alt(1:lx1,ix2,ix3),dim=1)
        if (x%alt(ix1apex,ix2,ix3)<0) error stop 'neutral background: field line is entirely below ground'
        if (lx1==1 .and. x%alt(1,ix2,ix3)>atmosbackground%altpmax) &
          error stop 'neutral background: insufficient altitude coverage'

        ! find the reference altitude and index beyond which the input data need to be extrapolated
        ! Start at the min(x1) end
        ix1=1
        do while (ix1<=ix1apex)
          if (x%alt(ix1,ix2,ix3)>atmosbackground%altpmax) exit
          ix1=ix1+1
        end do
        ix1ref1=ix1-1

        ! Now from the max(x1) end (in-case two-sided grid)
        ix1=lx1
        do while (ix1>=ix1apex)
          if (x%alt(ix1,ix2,ix3)>atmosbackground%altpmax) exit
          ix1=ix1-1
        end do
        ix1ref2=ix1+1

        !print*, ix1apex,ix1ref1,ix1ref2

        ! pull a reference density from the highest point on the simulation that exists below altitude limit
        !   of the input data
        ! Only extrapolate a hemisphere when the apex lies above file coverage.
        ! Two above-ground samples are required to determine its exponential slope.
        if (ix1apex>1 .and. ix1ref1<ix1apex) then
          if (ix1ref1<2) error stop 'neutral background: insufficient lower-end altitude coverage'
          if (x%alt(ix1ref1-1,ix2,ix3)<0) &
            error stop 'neutral background: insufficient above-ground lower-end samples'
          do ineu=1,lneu
            nrat=density_ratio(atmos%nnBG(ix1ref1,ix2,ix3,ineu),atmos%nnBG(ix1ref1-1,ix2,ix3,ineu))
            do ix1=ix1ref1+1,ix1apex
              atmos%nnBG(ix1,ix2,ix3,ineu)=nrat*atmos%nnBG(ix1-1,ix2,ix3,ineu)
            end do
          end do
        end if
        if (ix1apex<lx1 .and. ix1ref2>ix1apex) then
          if (ix1ref2>lx1-1) error stop 'neutral background: insufficient upper-end altitude coverage'
          if (x%alt(ix1ref2+1,ix2,ix3)<0) &
            error stop 'neutral background: insufficient above-ground upper-end samples'
          do ineu=1,lneu
            nrat=density_ratio(atmos%nnBG(ix1ref2,ix2,ix3,ineu),atmos%nnBG(ix1ref2+1,ix2,ix3,ineu))
            do ix1=ix1ref2-1,ix1apex,-1
              atmos%nnBG(ix1,ix2,ix3,ineu)=nrat*atmos%nnBG(ix1+1,ix2,ix3,ineu)
            end do
          end do
        end if

        ! temperature and velocity profiles get a ZOH
        if (ix1apex>1 .and. ix1ref1<ix1apex) then
          Tnref=atmos%TnBG(ix1ref1,ix2,ix3)
          vn1ref=atmos%vn1BG(ix1ref1,ix2,ix3)
          vn2ref=atmos%vn2BG(ix1ref1,ix2,ix3)
          vn3ref=atmos%vn3BG(ix1ref1,ix2,ix3)
          do ix1=ix1ref1+1,ix1apex
            atmos%TnBG(ix1,ix2,ix3)=Tnref
            atmos%vn1BG(ix1,ix2,ix3)=vn1ref
            atmos%vn2BG(ix1,ix2,ix3)=vn2ref
            atmos%vn3BG(ix1,ix2,ix3)=vn3ref
          end do
        end if
        if (ix1apex<lx1 .and. ix1ref2>ix1apex) then
          Tnref=atmos%TnBG(ix1ref2,ix2,ix3)
          vn1ref=atmos%vn1BG(ix1ref2,ix2,ix3)
          vn2ref=atmos%vn2BG(ix1ref2,ix2,ix3)
          vn3ref=atmos%vn3BG(ix1ref2,ix2,ix3)
          do ix1=ix1ref2-1,ix1apex,-1
            atmos%TnBG(ix1,ix2,ix3)=Tnref
            atmos%vn1BG(ix1,ix2,ix3)=vn1ref
            atmos%vn2BG(ix1,ix2,ix3)=vn2ref
            atmos%vn3BG(ix1,ix2,ix3)=vn3ref
          end do
        end if

        ! any addition 'fixing' that needs to be done on this field-line profile...
        atmos%vn1BG(1:lx1,ix2,ix3)=atmos%vn1BG(1:lx1,ix2,ix3)*(0.5 + 0.5*tanh((x%alt(1:lx1,ix2,ix3)-150e3)/10e3))

        !! Now we need to find the below ground points and do something with them
        ix1=1
        do while (ix1<=ix1apex)
          if (x%alt(ix1,ix2,ix3)>=0) exit
          ix1=ix1+1
        end do
        ix1ref1=max(ix1-1,1)

        ix1=lx1
        do while (ix1>=ix1apex)
          if (x%alt(ix1,ix2,ix3)>=0) exit
          ix1=ix1-1
        end do
        ix1ref2=min(ix1+1,lx1)

        !print*,ix1ref1,ix1ref2

        do ineu=1,lneu
          if (x%alt(1,ix2,ix3)<0) then
            do ix1=1,ix1ref1
              atmos%nnBG(ix1,ix2,ix3,ineu)=atmos%nnBG(ix1ref1+1,ix2,ix3,ineu)
            end do
          end if
          if (x%alt(lx1,ix2,ix3)<0) then
            do ix1=lx1,ix1ref2,-1
              atmos%nnBG(ix1,ix2,ix3,ineu)=atmos%nnBG(ix1ref2-1,ix2,ix3,ineu)
            end do
          end if
        end do

        if (x%alt(1,ix2,ix3)<0) then
          do ix1=1,ix1ref1
            atmos%TnBG(ix1,ix2,ix3)=atmos%TnBG(ix1ref1+1,ix2,ix3)
          end do
        end if
        if (x%alt(lx1,ix2,ix3)<0) then
          do ix1=lx1,ix1ref2,-1
            atmos%TnBG(ix1,ix2,ix3)=atmos%TnBG(ix1ref2-1,ix2,ix3)
          end do
        end if
      end do
    end do

  end subroutine neutral_background_fileinput_copyout


  real(wp) function density_ratio(nref,nlower) result(ratio)
    real(wp), intent(in) :: nref,nlower

    if (.not. ieee_is_finite(nref) .or. .not. ieee_is_finite(nlower)) &
      error stop 'neutral background: nonfinite extrapolation density'
    if (nref<0 .or. nlower<0) error stop 'neutral background: negative extrapolation density'
    ! A species absent at the upper sample remains absent above it. A positive
    ! upper density over a zero lower density cannot define an exponential.
    if (nref==0) then
      ratio=0
    else
      if (nlower==0) error stop 'neutral background: zero lower extrapolation density'
      ratio=nref/nlower
    end if
  end function density_ratio


  !> initializes neutral atmosphere by:
  !    1)  establishing initial background for density, temperature, and winds from MSIS and HWM
  !! Arguably this should be called for init and for neutral updates, except for the allocation part...
  subroutine neutral_background_empirical(cfg,ymd,UTsec,x,v2grid,v3grid,atmos)
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
    !call cpu_time(tstart)
    call neutral_atmos(ymd,UTsec,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3), &
                         cfg%activ,cfg%msis_version,atmos)
    !if (mpi_cfg%myid == 0) then
    !  call cpu_time(tfin)
      !print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if

    !> Horizontal wind model initialization/background
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    !call cpu_time(tstart)
    call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x, atmos=atmos)
    !! we sum the horizontal wind with the background state vector
    !! if HWM14 is disabled, neutral_winds returns the background state vector unmodified
    !if (mpi_cfg%myid == 0) then
    !  call cpu_time(tfin)
      !print *, 'Initial neutral winds (from HWM) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if

    !call neutral_aggregate(v2grid,v3grid,atmos,atmosperturb,.true.)    !.true. gives just background state at start
  end subroutine neutral_background_empirical


  !> initialization procedure needed for MSIS 2.0
  subroutine msisinit_in(cfg)
    type(gemini_cfg), intent(in) :: cfg
    logical :: exists

    character(len=11) :: msis2_param_file

    select case (cfg%msis_version)
    case(0)
      !! MSISE00
      return
    case(20)
      msis2_param_file = "msis20.parm"
    case(21)
      msis2_param_file = "msis21.parm"
    case default
      !! new or unknown version of MSIS, default MSIS 2.x parameter file
      msis2_param_file = ""
    end select

    if(len_trim(msis2_param_file) > 0) then
      inquire(file=msis2_param_file, exist=exists)
      if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
        ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
        'This limitation comes from how MSIS 2.x is coded internally.'
      call msisinit(parmfile=msis2_param_file)
    else
      call msisinit()
    end if
  end subroutine msisinit_in


  subroutine clear_neutral_background_fileinput(atmosbackground)
    type(neutraldataBG), pointer, intent(inout) :: atmosbackground

    if (associated(atmosbackground)) deallocate(atmosbackground)    ! should nullify atmosbackground

  end subroutine clear_neutral_background_fileinput
end module neutral_background
