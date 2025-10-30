module neutraldataBGobj

! type extension for file-based neutral background atmospheric data from a profile

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,debug,pi
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize1,get_grid1!,get_neutralBG
use timeutils, only: dateinc,date_filename

implicit none (type, external)
private
public :: neutralBGdata

type, extends(inputdata) :: neutralBGdata
  ! coordinates for input precipitation data, and storage
  real(wp), dimension(:), pointer :: altp
  integer, pointer :: llon,llat,lalt
  real(wp), dimension(:,:), pointer :: natmp
  real(wp), dimension(:,:,:,:), pointer :: natmiprev,natminext,natminow

  ! work and target coordinates
  real(wp), dimension(:,:,:), allocatable :: altimat
  real(wp), dimension(:), pointer :: alti

  contains
    ! deferred bindings
    procedure :: init=>init_neutralBG
    procedure :: set_coordsi=>set_coordsi_neutralBG
    procedure :: load_data=>load_data_neutralBG
    procedure :: load_grid=>load_grid_neutralBG
    procedure :: load_size=>load_size_neutralBG     ! load the size of the input data files
    final :: destructor
end type neutralBGdata

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_neutralBG(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutralBGdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character(*), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtmodel,dtdata                      ! model time step and cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname

    ! need to allow interpolation from 2D to 3D
    !self%flagallow2D3D=.true.

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='neutral background'
    call self%set_name(strname)
    self%flagdoinput=cfg%flagneutralBGfile/=0

    ! read the simulation size from the source directory and allocate arrays
    allocate(self%lc1,self%lc2,self%lc3)      ! these are pointers
    self%llon=>self%lc2;  self%llat=>self%lc3;  self%lalt=>self%lc1;
    call self%load_size()
    call self%set_sizes(0, &
                       0,0,0, &
                       0,0,0, &
                       6, &     ! target data for neutralBG info is a 3D set of arrays
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%altp=>self%coord1
    call self%load_grid()

    ! Set input data array pointers to faciliate easy to read input code; these may or may not be helpful to user
    self%natmp=>self%data3D(:,1,1,:)   !contiguous in memory since only one element in x2,3
    self%natmiprev=>self%data3Di(:,:,:,:,1)     
    self%natminext=>self%data3Di(:,:,:,:,2)
    self%natminow=>self%data3Dinow(:,:,:,:)

    ! must initialize prev state or else the first set of data will not be interpolated correctly
    self%natmiprev=0.0
    self%natminext=0.0

    ! set to start time of simulation - will be set first time update is called
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_neutralBG


  !> get the input grid size from file, all workers will just call this sicne this is a one-time thing
  subroutine load_size_neutralBG(self)
    class(neutralBGdata), intent(inout) :: self
    integer :: ltmp    ! throwaway variable

    ! basic error checking
    if (.not. self%flagsource) error stop 'neutralBGdata:load_size_neutralBG() - must define a source directory first'

    ! read sizes
    !print '(/,A,/,A)', 'neutralBG input:','--------------------'
    !print '(A)', 'READ neutralBG size from: ' // self%sourcedir
    call get_simsize1(self%sourcedir // "/simsize.h5", self%lalt)   ! mangling lat->alt

    self%llon=1; self%llat=1;    ! force input to be only profile-based

    !print '(A,2I6)', 'soflux size: llon,llat:  ',self%llon,self%llat
    if (self%lalt < 1) then
     print*, '  neutralBG grid size must be strictly positive: ' //  self%sourcedir
     error stop
    end if

    ! flag to denote input data size is set
    self%flagdatasize=.true.
  end subroutine load_size_neutralBG


  !> get the grid information from a file, all workers will just call this since one-time
  subroutine load_grid_neutralBG(self)
    class(neutralBGdata), intent(inout) :: self

    ! read grid data
    call get_grid1(self%sourcedir // "/simgrid.h5", self%altp)

    !print '(A,4F9.3)', 'Solar flux glon,glat extent:  ',minval(self%glonp(:)),maxval(self%glonp(:)), &
    !                                                       minval(self%glatp(:)),maxval(self%glatp(:))
    if(.not. all(ieee_is_finite(self%altp))) error stop 'neutralBGBCs_fileinput: glon must be finite'
    !if(.not. all(ieee_is_finite(self%glatp))) error stop 'neutralBGBCs_fileinput: glat must be finite'

    !print*, 'min/max glonp:  ',minval(self%glonp),maxval(self%glonp)
    !print*, 'min/max glatp:  ',minval(self%glatp),maxval(self%glatp)
  end subroutine load_grid_neutralBG


  !> set target coordinates for interpolation sights
  subroutine set_coordsi_neutralBG(self,cfg,x)
    class(neutralBGdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg     ! presently not used but possibly eventually?
    class(curvmesh), intent(in) :: x
    integer :: ix1,ix2,ix3

    ! aliases for target interpolation sites
    self%alti=>self%coord1i

    allocate(self%altimat(1:x%lx1,1:x%lx2,1:x%lx3))     ! why not local variables?  FIXME

    ! Target coordinates are 3D in this case, e.g. due to dipole grid where alt spacing varies across domain
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          self%altimat(ix1,ix2,ix3)=x%alt(ix1,ix2,ix3)
        end do
      end do
    end do
    self%alti=pack(self%altimat,.true.)

    deallocate(self%altimat)
    self%flagcoordsi=.true.

    !print*, 'min/max gloni:  ',minval(self%gloni),maxval(self%gloni)
    !print*, 'min/max glati:  ',minval(self%glati),maxval(self%glati)
  end subroutine set_coordsi_neutralBG


  !> have all processes read in data from file to avoid any message passing
  subroutine load_data_neutralBG(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutralBGdata), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warnings

    !! all workers should update the date
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt, ymdtmp, UTsectmp)

    !!!!!!  read in solar flux data from file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! this read must be done repeatedly through simulation so have only root do file io
    !print*, '  date and time:  ',ymdtmp,UTsectmp
    !print*, '  neutralBG filename:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp)
    ! read in the data for the "next" frame from file
    !call get_neutralBG(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", self%Iinfp)
    !print*, 'min/max data:  ',  minval(self%Iinfp),maxval(self%Iinfp)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !print*, ' precip data succesfully input...'
  end subroutine load_data_neutralBG


  !> destructor needs to clear memory out
  subroutine destructor(self)
    type(neutralBGdata), intent(inout) :: self

    call self%dissociate_pointers()
  end subroutine destructor
end module neutraldataBGobj
