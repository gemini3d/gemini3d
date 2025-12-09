module neutral_background

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
        clear_neutral_background_fileinput

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
    logical :: flagnoninverted, flagtwosided
    real(wp) :: nref,Tnref,vn1ref,vn2ref,vn3ref
    real(wp) :: altimax1,altimax2,altapex,nrat

    ! place the background atmospheric density and temperature (drifts still need to be rotated)
    !   from object into arrays used by other modules
    atmos%nnBG(:,:,:,1:5)=atmosbackground%natminow(:,:,:,1:5)
    atmos%TnBG=atmosbackground%natminow(:,:,:,9)
    call NO_calc(atmos%nnBG,atmos%TnBG)
    
    ! rotate winds to model native and place into wind background arrays in the atmos object
    call rotate_geo2native(atmosbackground%natminow(:,:,:,6), &
            atmosbackground%natminow(:,:,:,7), &
            atmosbackground%natminow(:,:,:,8), &
            x,atmos,.true.)

    flagnoninverted=x%alt(2,1,1)>x%alt(1,1,1)
    flagtwosided=abs(x%alt(1,1,1)-x%alt(lx1,1,1))<1.0e3_wp
    lneu=size(atmos%nnBG,4)

    ! FIXME: assumes monotonic altitude array and overwrites if not!!!
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        ! for each field line we need to find the apex altitude and index (it possibly could be a different index
        !   on each field line???
        altapex=maxval(x%alt(:,ix2,ix3))
        ix1apex=maxloc(x%alt(:,ix2,ix3),dim=1)

        ! find the reference altitude and index beyond which the input data need to be extrapolated
        ! Start at the min(x1) end
        altimax1=0.0
        ix1ref1=1
        do ix1=1,ix1apex    ! we only iterate up to what we know is the apex altitude for this field line
          if (x%alt(ix1,ix2,ix3)<atmosbackground%altpmax .and. x%alt(ix1,ix2,ix3)>altimax1) then
            altimax1=x%alt(ix1,ix2,ix3)
            ix1ref1=ix1
          end if
        end do

        ! Now from the max(x1) end (in-case two-sided grid)
        altimax2=0.0
        ix1ref2=lx1
        do ix1=lx1,ix1apex,-1
          if (x%alt(ix1,ix2,ix3)<atmosbackground%altpmax .and. x%alt(ix1,ix2,ix3)>altimax2) then
            altimax2=x%alt(ix1,ix2,ix3)
            ix1ref2=ix1
          end if         
        end do

        ! pull a reference density from the highest point on the simulation that exists below altitude limit 
        !   of the input data
        do ineu=1,lneu
          nref=atmos%nnBG(ix1ref1,ix2,ix3,ineu)
          !if (flagnoninverted) then
            nrat=nref/atmos%nnBG(ix1ref1-1,ix2,ix3,ineu)
            do ix1=ix1ref1+1,ix1apex
              if (ix1>1) then
                atmos%nnBG(ix1,ix2,ix3,ineu)=nrat*atmos%nnBG(ix1-1,ix2,ix3,ineu)
              else
                error stop 'problem extrapolating background neutral density profile:  index OOB'
              end if
            end do
          !else

          nref=atmos%nnBG(ix1ref2,ix2,ix3,ineu)
            nrat=nref/atmos%nnBG(ix1ref2+1,ix2,ix3,ineu)
            do ix1=ix1ref2-1,ix1apex,-1
              if (ix1<lx1) then
                atmos%nnBG(ix1,ix2,ix3,ineu)=nrat*atmos%nnBG(ix1+1,ix2,ix3,ineu)
              else
                error stop 'problem extrapolating background neutral density profile.  index OOB'               
              end if
            end do
          !end if
        end do

        ! temperature and velocity profiles get a ZOH
        Tnref=atmos%TnBG(ix1ref1,ix2,ix3)
        vn1ref=atmos%vn1(ix1ref1,ix2,ix3)
        vn2ref=atmos%vn2(ix1ref1,ix2,ix3)
        vn3ref=atmos%vn3(ix1ref1,ix2,ix3)
        !if (flagnoninverted) then
          do ix1=ix1ref1+1,ix1apex
            atmos%TnBG(ix1,ix2,ix3)=Tnref
            atmos%vn1BG(ix1,ix2,ix3)=vn1ref
            atmos%vn2BG(ix1,ix2,ix3)=vn2ref
            atmos%vn3BG(ix1,ix2,ix3)=vn3ref
          end do
        !else
        Tnref=atmos%TnBG(ix1ref2,ix2,ix3)
        vn1ref=atmos%vn1(ix1ref2,ix2,ix3)
        vn2ref=atmos%vn2(ix1ref2,ix2,ix3)
        vn3ref=atmos%vn3(ix1ref2,ix2,ix3)
          do ix1=ix1ref2-1,ix1apex,-1
            atmos%TnBG(ix1,ix2,ix3)=Tnref
            atmos%vn1BG(ix1,ix2,ix3)=vn1ref
            atmos%vn2BG(ix1,ix2,ix3)=vn2ref
            atmos%vn3BG(ix1,ix2,ix3)=vn3ref
          end do
        !end if

        ! any addition 'fixing' that needs to be done on this field-line profile...
        atmos%vn1BG(1:lx1,ix2,ix3)=atmos%vn1BG(1:lx1,ix2,ix3)*(0.5 + 0.5*tanh((x%alt(1:lx1,ix2,ix3)-150e3)/10e3))
!        atmos%vn2BG(1:lx1,ix2,ix3)=atmos%vn2BG(1:lx1,ix2,ix3)*(0.5 + 0.5*tanh((x%alt(1:lx1,ix2,ix3)-150e3)/10e3))
!        atmos%vn3BG(1:lx1,ix2,ix3)=atmos%vn3BG(1:lx1,ix2,ix3)*(0.5 + 0.5*tanh((x%alt(1:lx1,ix2,ix3)-150e3)/10e3))

        ix1ref1=1
        do ix1=1,ix1apex    ! we only iterate up to what we know is the apex altitude for this field line
          if (x%alt(ix1,ix2,ix3)<0) then
            ix1ref1=ix1
          end if
        end do
        ix1ref2=lx1
        do ix1=lx1,ix1apex,-1
          if (x%alt(ix1,ix2,ix3)<0) then
            ix1ref2=ix1
          end if         
        end do

        do ineu=1,lneu
          do ix1=1,ix1ref1
            atmos%nnBG(ix1,ix2,ix3,ineu)=atmos%nnBG(ix1ref1+1,ix2,ix3,ineu)
          end do
          do ix1=lx1,ix1ref2,-1
            atmos%nnBG(ix1,ix2,ix3,ineu)=atmos%nnBG(ix1ref2-1,ix2,ix3,ineu)           
          end do
        end do

        do ix1=1,ix1ref1
          atmos%TnBG(ix1,ix2,ix3)=atmos%TnBG(ix1ref1+1,ix2,ix3)
        end do
        do ix1=lx1,ix1ref2,-1
          atmos%TnBG(ix1,ix2,ix3)=atmos%TnBG(ix1ref2-1,ix2,ix3)           
        end do
      end do
    end do

  end subroutine neutral_background_fileinput_copyout


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
    call neutral_atmos(cfg%ymd0,cfg%UTsec0,x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3), &
                         cfg%activ,cfg%msis_version,atmos)
    !if (mpi_cfg%myid == 0) then
    !  call cpu_time(tfin)
      !print *, 'Initial neutral density and temperature (from MSIS) at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
    !end if

    !> Horizontal wind model initialization/background
    !if (mpi_cfg%myid == 0) call cpu_time(tstart)
    !call cpu_time(tstart)
    call neutral_winds(cfg%ymd0, cfg%UTsec0, Ap=cfg%activ(3), x=x, atmos=atmos)
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
