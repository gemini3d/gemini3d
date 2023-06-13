module neutraldata2Dobj

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_grid2,get_neutral2
use timeutils, only: dateinc,date_filename
use h5fortran, only: hdf5_file
use grid, only: gridflag

implicit none (type,external)
private
public :: neutraldata2D

!> type definition for 3D neutral data
type, extends(neutraldata), abstract :: neutraldata2D
  ! source data coordinate pointers
  real(wp), dimension(:), pointer :: horzn,zn
  integer, pointer :: lhorzn,lxn,lzn        ! lxn not used but still bound

  ! work arrays needed by various procedures re: target coordinates
  real(wp), dimension(:,:,:), allocatable :: horzimat,zimat
  real(wp), dimension(:), pointer :: zi,horzi

  ! source data pointers, note we only have *horizontal* wind components here
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvnz,dvnhorz,dTn

  ! projection factors needed to rotate input data onto grid
  real(wp), dimension(:,:,:), allocatable :: proj_ezp_e1,proj_ezp_e2,proj_ezp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_ehorzp_e1,proj_ehorzp_e2,proj_ehorzp_e3
  contains
    ! replacement for gridsize and gridload
    procedure(loadneu2D), deferred :: load_sizeandgrid_neu2D
    procedure :: rotate_winds
    procedure :: init_neu2D_simple

    ! overriding procedures
    procedure :: update
    procedure :: init_storage

    ! bindings for deferred procedures
    procedure :: load_data=>load_data_neu2D
    procedure :: load_grid=>load_grid_neu2D    ! stub, does nothing see load_sizeandgrid_neu3D()
    procedure :: load_size=>load_size_neu2D    ! stub, does nothing "

    ! defer to type extensions
    !procedure :: set_coordsi=>set_coordsi_neu2D

    ! clean up some of the pointers
    procedure :: dissociate_neutral2D_pointers
end type neutraldata2D

interface
  subroutine loadneu2D(self,cfg)
    import neutraldata2D,gemini_cfg
    class(neutraldata2D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
  end subroutine loadneu2D
end interface

contains
  !> initialize storage for this type of neutral input data
  subroutine init_neu2D_simple(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutraldata2D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel,dtdata
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname    ! allow auto-allocate for strings

    ! force 3D interpolation regardless of working subarray size
    !self%flagforcenative=.true.

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='neutral perturbations (2D)'
    call self%set_name(strname)
    call self%set_cadence(dtdata)
    self%flagdoinput=cfg%flagdneu/=0

    ! set sizes, we have 7 arrays all 3D (irrespective of 2D vs. 3D neutral input).  for 3D neutral input
    !    the situation is more complicated that for other datasets because you cannot compute the number of
    !    source grid points for each worker until you have root compute the entire grid and slice everything up
    allocate(self%lc1,self%lc2,self%lc3)                                     ! these are pointers, even though scalar
    !! this is a bit tricky because the inpudata class wants non-singleton dimension to be the same; here lc1,lc2
    self%lzn=>self%lc1; self%lhorzn=>self%lc2;
    self%lxn=>self%lc3                             ! these referenced while reading size and grid data
    call self%set_coordsi(cfg,x)                   ! since this preceeds init_storage it must do the work of allocating some spaces
    call self%load_sizeandgrid_neu2D(cfg)          ! cfg needed to form source neutral grid
    call self%set_sizes( &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data along pairs of axes
             6, &          ! number 3D datasets, for neutraldata2D we have singleton dimensions for 2D input
             x)            ! The main purpose of this is to set the number of 3D datasets (other params already set)

    ! allocate space for arrays, note for neutrals some of this has already happened so there is an overloaded procedure
    call self%init_storage()

    ! set aliases to point to correct source data arrays
    self%dnO=>self%data3D(:,:,:,1)
    self%dnN2=>self%data3D(:,:,:,2)
    self%dnO2=>self%data3D(:,:,:,3)
    self%dvnz=>self%data3D(:,:,:,4)
    self%dvnhorz=>self%data3D(:,:,:,5)
    self%dTn=>self%data3D(:,:,:,6)    ! note this gets interpolated into a variable bound to vn3i, so need to adjust during rotations

    ! call to base class procedure to set pointers for prev,now,next (must already have space allocated)
    call self%setptrs_grid()

    ! initialize previous data so we get a correct starting value
    self%dnOiprev=0
    self%dnN2iprev=0
    self%dnO2iprev=0
    self%dvn1iprev=0
    self%dvn2iprev=0
    self%dvn3iprev=0
    self%dTniprev=0

    ! set to start time of simulation - not needed since assigned by update on first call.  FIXME: a bit messy
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_neu2D_simple


  !> create storage for arrays needed specifically for 3D neutral input calculations, overrides the base class procedure
  subroutine init_storage(self)
    class(neutraldata2D), intent(inout) :: self
    integer :: lc1,lc2,lc3
    integer :: lc1i,lc2i,lc3i
    integer :: l0D
    integer :: l1Dax1,l1Dax2,l1Dax3
    integer :: l2Dax23,l2Dax12,l2Dax13
    integer :: l3D

    ! check sizes are set
    if (.not. self%flagsizes) error stop 'inpudata:init_storage(); must set sizes before allocations...'

    ! local size variables for convenience
    lc1=self%lc1; lc2=self%lc2; lc3=self%lc3;
    lc1i=self%lc1i; lc2i=self%lc2i; lc3i=self%lc3i;
    l0D=self%l0D
    l1Dax1=self%l1Dax1; l1Dax2=self%l1Dax2; l1Dax3=self%l1Dax3;
    l2Dax23=self%l2Dax23; l2Dax12=self%l2Dax12; l2Dax13=self%l2Dax13;
    l3D=self%l3D

    ! NOTE: type extensions are reponsible for zeroing out any arrays they will use in their own init() bindings

    ! input data coordinate arrays are set by load_gridandsize()

    ! allocate target coords, for neutral3D the standard set (coord1i, etc.) are done in set_coordsi()
    allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
    allocate(self%coord2iax23(lc2i*lc3i),self%coord3iax23(lc2i*lc3i))
    allocate(self%coord1iax13(lc1i*lc3i),self%coord3iax13(lc1i*lc3i))
    allocate(self%coord1iax12(lc1i*lc2i),self%coord2iax12(lc1i*lc2i))

    ! allocate object arrays for input data at a reference time.  FIXME: do we even need to store this perm. or can be local to
    ! load_data?
    allocate(self%data0D(l0D))
    allocate(self%data1Dax1(lc1,l1Dax1), self%data1Dax2(lc2,l1Dax2), self%data1Dax3(lc3,l1Dax3))
    allocate(self%data2Dax23(lc2,lc3,l2Dax23), self%data2Dax12(lc1,lc2,l2Dax12), self%data2Dax13(lc1,lc3,l2Dax13))
    allocate(self%data3D(lc1,lc2,lc3,l3D))

    ! allocate object arrays for interpolation sites at reference times
    allocate(self%data0Di(l0D,2))
    allocate(self%data1Dax1i(lc1i,l1Dax1,2), self%data1Dax2i(lc2i,l1Dax2,2), self%data1Dax3i(lc3i,l1Dax3,2))
    allocate(self%data2Dax23i(lc2i,lc3i,l2Dax23,2), self%data2Dax12i(lc1i,lc2i,l2Dax12,2), self%data2Dax13i(lc1i,lc3i,l2Dax13,2))
    allocate(self%data3Di(lc1i,lc2i,lc3i,l3D+1,2))     ! note extra parameter for three vector components!!!

    ! allocate object arrays at interpolation sites for current time.  FIXME: do we even need to store permanently?
    allocate(self%data0Dinow(l0D))
    allocate(self%data1Dax1inow(lc1i,l1Dax1), self%data1Dax2inow(lc2i,l1Dax2), self%data1Dax3inow(lc3i,l1Dax3))
    allocate(self%data2Dax23inow(lc2i,lc3i,l2Dax23), self%data2Dax12inow(lc1i,lc2i,l2Dax12), self%data2Dax13inow(lc1i,lc3i,l2Dax13))
    allocate(self%data3Dinow(lc1i,lc2i,lc3i,l3D+1))    ! +1 because even with 2D input we still need to track 3 comps.

    self%flagalloc=.true.
  end subroutine init_storage


  !> do nothing stub
  subroutine load_size_neu2D(self)
    class(neutraldata2D), intent(inout) :: self

    integer :: i
    i = self%lxn
    !! avoid unused argument warnings

  end subroutine load_size_neu2D


  !> do nothing stub
  subroutine load_grid_neu2D(self)
    class(neutraldata2D), intent(inout) :: self

    integer :: i
    i = self%lxn
    !! avoid unused argument warnings

  end subroutine load_grid_neu2D


  !> Have all workers separately load data out of file to avoid message passing
  subroutine load_data_neu2D(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutraldata2D), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
    integer :: lhorzn,lzn                        !number of horizontal grid points

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warnings

    ! sizes for convenience
    lhorzn=self%lhorzn; lzn=self%lzn;

    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt,ymdtmp,UTsectmp)                !get the date for "next" params

    call get_neutral2(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", &
      self%dnO,self%dnN2,self%dnO2,self%dvnhorz,self%dvnz,self%dTn)

    !print*, 'Loading 2D neutral data from:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5"

    if (debug) then
      print *, 'Min/max values for dnO:  ',minval(self%dnO),maxval(self%dnO)
      print *, 'Min/max values for dnN:  ',minval(self%dnN2),maxval(self%dnN2)
      print *, 'Min/max values for dnO:  ',minval(self%dnO2),maxval(self%dnO2)
      print *, 'Min/max values for dvnhorz:  ',minval(self%dvnhorz),maxval(self%dvnhorz)
      print *, 'Min/max values for dvnz:  ',minval(self%dvnz),maxval(self%dvnz)
      print *, 'Min/max values for dTn:  ',minval(self%dTn),maxval(self%dTn)
    endif

    if (.not. all(ieee_is_finite(self%dnO))) error stop 'dnO: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%dnN2))) error stop 'dnN2: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%dnO2))) error stop 'dnO2: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%dvnhorz))) error stop 'dvnhorz: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%dvnz))) error stop 'dvnz: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%dTn))) error stop 'dTn: non-finite value(s)'

    ! print some diagnostics for the input data
    if (debug) then
      print*, 'neutral data size:  ',lhorzn,lzn
      print *, 'Min/max values for dnO:  ',minval(self%dnO),maxval(self%dnO)
      print *, 'Min/max values for dnN:  ',minval(self%dnN2),maxval(self%dnN2)
      print *, 'Min/max values for dnO:  ',minval(self%dnO2),maxval(self%dnO2)
      print *, 'Min/max values for dvnrho:  ',minval(self%dvnhorz),maxval(self%dvnhorz)
      print *, 'Min/max values for dvnz:  ',minval(self%dvnz),maxval(self%dvnz)
      print *, 'Min/max values for dTn:  ',minval(self%dTn),maxval(self%dTn)
      !print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
    end if
  end subroutine load_data_neu2D


  !> overriding procedure for updating neutral atmos (need additional rotation steps)
  subroutine update(self,cfg,dtmodel,t,x,ymd,UTsec)
    class(neutraldata2D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel             ! need both model and input data time stepping
    real(wp), intent(in) :: t                   ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x            ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations

    ! execute a basic update
    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)

    ! FIXME: more efficient to rotate only when a read/interpolate in space is done?
    ! now we need to rotate velocity fields following interpolation (they are magnetic ENU prior to this step)
    call self%rotate_winds()

    ! print some diagnostic data once the udpate has occurred
    if (debug) then
      print*, ''
      print*, 'neutral data size:  ',self%lzn,self%lhorzn,self%lxn
      print*, 'neutral data time:  ',ymd,UTsec
      print*, ''
      print *, 'Min/max values for dnOinext:  ',minval(self%dnOinext),maxval(self%dnOinext)
      print *, 'Min/max values for dnN2inext:  ',minval(self%dnN2inext),maxval(self%dnN2inext)
      print *, 'Min/max values for dnO2inext:  ',minval(self%dnO2inext),maxval(self%dnO2inext)
      print *, 'Min/max values for dvn1inext:  ',minval(self%dvn1inext),maxval(self%dvn1inext)
      print *, 'Min/max values for dvn2inext:  ',minval(self%dvn2inext),maxval(self%dvn2inext)
      print *, 'Min/max values for dvn3inext:  ',minval(self%dvn3inext),maxval(self%dvn3inext)
      print *, 'Min/max values for dTninext:  ',minval(self%dTninext),maxval(self%dTninext)
      print*, ''
      print *, 'Min/max values for dnOinow:  ',minval(self%dnOinow),maxval(self%dnOinow)
      print *, 'Min/max values for dnN2inow:  ',minval(self%dnN2inow),maxval(self%dnN2inow)
      print *, 'Min/max values for dnO2inow:  ',minval(self%dnO2inow),maxval(self%dnO2inow)
      print *, 'Min/max values for dvn1inow:  ',minval(self%dvn1inow),maxval(self%dvn1inow)
      print *, 'Min/max values for dvn2inow:  ',minval(self%dvn2inow),maxval(self%dvn2inow)
      print *, 'Min/max values for dvn3inow:  ',minval(self%dvn3inow),maxval(self%dvn3inow)
      print *, 'Min/max values for dTninow:  ',minval(self%dTninow),maxval(self%dTninow)
    end if
  end subroutine update


  !> This subroutine takes winds stored in self%dvn?inow and applies a rotational transformation onto the
  !      grid object for this simulation.  Provided that the horizontal projections have been computed
  !      correctly the same rotation can be used for axisymmetric and cartesian.
  subroutine rotate_winds(self)
    class(neutraldata2D), intent(inout) :: self
    integer :: ix1,ix2,ix3
    real(wp) :: vnhorz,vnz,Tn

    ! do rotations one grid point at a time to cut down on temp storage needed.  Note that until this point there
    !   shoudl be only zero data stored in vn3 since this class is for 2D data input, instead temperature
    !   gets stored in the dvn3i variables.
    do ix3=1,self%lc3i
      do ix2=1,self%lc2i
        do ix1=1,self%lc1i
          vnz=self%dvn1inow(ix1,ix2,ix3)
          vnhorz=self%dvn2inow(ix1,ix2,ix3)
          Tn=self%dvn3inow(ix1,ix2,ix3)    ! need to save because it will get overwritten in rotation
          self%dvn1inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e1(ix1,ix2,ix3) + &
                                        vnhorz*self%proj_ehorzp_e1(ix1,ix2,ix3)
          self%dvn2inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e2(ix1,ix2,ix3) + &
                                        vnhorz*self%proj_ehorzp_e2(ix1,ix2,ix3)
          self%dvn3inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e3(ix1,ix2,ix3) + &
                                        vnhorz*self%proj_ehorzp_e3(ix1,ix2,ix3)
          self%dTninow(ix1,ix2,ix3)=Tn     ! assign saved temperature into correct slot in "output" variables
        end do
      end do
    end do
  end subroutine rotate_winds

  !> destructor for when object goes out of scope
  subroutine dissociate_neutral2D_pointers(self)
    class(neutraldata2D) :: self

    ! deallocate arrays from base inputdata class
    call self%dissociate_pointers()

    ! null pointers specific to parent neutraldata class
    call self%dissociate_neutral_pointers()

    ! now deallocate arrays specific to this extension
    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
    deallocate(self%proj_ehorzp_e1,self%proj_ehorzp_e2,self%proj_ehorzp_e3)
    deallocate(self%horzimat,self%zimat)

    ! set pointers to null
    nullify(self%horzi,self%zi);
    nullify(self%horzn,self%zn);
    nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnhorz,self%dTn)
  end subroutine dissociate_neutral2D_pointers
end module neutraldata2Dobj
