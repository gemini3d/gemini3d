module inputdataobj

use phys_consts, only : wp
use gemini3d_config, only: gemini_cfg
use meshobj, only : curvmesh
use meshobj_dipole, only: dipolemesh
use interpolation, only : interp1,interp2,interp3
use timeutils, only : dateinc, date_filename, find_lastdate

implicit none (type, external)
private
public :: inputdata


!> this is a generic class for an data object being input into the model and interpolated in space and time
type, abstract :: inputdata
  character(:), allocatable :: dataname     ! string description of dataset
  character(:), allocatable :: sourcedir    ! source location containing data input files

  !! flags to aid error checking
  logical :: flagdatasize=.false.       ! input data sizes set?
  logical :: flagsizes=.false.          ! all sizes set?
  logical :: flagalloc=.false.          ! space for data allocated?
  logical :: flagprimed=.false.         ! initial setup of input data files (priming)?
  logical :: flagcadence=.false.        ! time cadence of input has been set?
  logical :: flagsource=.false.         ! source directory for data set?
  logical :: flagcoordsi=.false.        ! interpolation sites set?
  logical :: flagforcenative=.false.    ! force all interpolations to be done with native array rank rather than detecting singleton
  logical :: flagallow2D3D=.false.      ! allow a dataset to interpolate 2D to 3D
  logical :: flagdoinput=.false.        ! extensions need to define how they know whether or not they need to do file input
  logical :: flagfirst=.true.           ! true prior to performing first update
  logical :: flagdipmesh=.false.        ! are we interpolating to a dipole mesh?

  !! here we store data that have already been received but not yet interpolated
  real(wp), dimension(:), pointer :: coord1,coord2,coord3     ! coordinates for the source data (interpolant coords)
  integer, pointer :: lc1,lc2,lc3                                      ! dataset length along the 3 coordinate axes
  real(wp), dimension(:), pointer :: data0D
  real(wp), dimension(:,:), pointer :: data1Dax1,data1Dax2,data1Dax3
  real(wp), dimension(:,:,:), pointer :: data2Dax23,data2Dax12,data2Dax13
  real(wp), dimension(:,:,:,:), pointer :: data3D

  !! here we store data that have already been spatially interpolated
  real(wp), dimension(:,:), pointer :: data0Di                    ! array for storing a "stack" of scalar data (only interpolated in time)
                                                                     !  last axis is for prev,next copies of data for interp in time
                                                                     !  second to last axis is for number of datasets of this dimension
  integer :: l0D                                                     ! length/number of scalar datasets
  real(wp), dimension(:,:,:), pointer :: data1Dax1i                  ! array for storing series of 1D data, array varies along non-singleton axis
  real(wp), dimension(:,:,:), pointer :: data1Dax2i,data1Dax3i       ! 1D data arrays varying along coordinates (axes) 2 and 3
  integer :: l1Dax1,l1Dax2,l1Dax3
  real(wp), dimension(:,:,:,:), pointer :: data2Dax23i               ! array for storing series of 2D data, varies along two non-singleton axes
  real(wp), dimension(:,:,:,:), pointer :: data2Dax12i,data2dax13i   !2D arrays varying along 1,2 and 1,3 axes
  integer :: l2Dax23,l2Dax12,l2Dax13
  real(wp), dimension(:,:,:,:,:), pointer :: data3Di                 ! array for storing series of 3D data
  integer :: l3D

  !! by default we have one set of target coordinates; extension can define others, if needed; note these are "flat" arrays (rank 1)
  real(wp), dimension(:), pointer :: coord1i,coord2i,coord3i             ! coordinates of the interpolation sites, full 3D, size lc1i*lc2i*lc3i
  real(wp), dimension(:), pointer :: coord1iax1                          ! 1D target coords for variations along axis 1
  real(wp), dimension(:), pointer :: coord2iax2                          ! 1D target coords for variations along axis 2
  real(wp), dimension(:), pointer :: coord3iax3                          ! 1D target coords for variations along axis 3
  real(wp), dimension(:), pointer :: coord2iax23,coord3iax23             ! 2D target along axes 2,3
  real(wp), dimension(:), pointer :: coord1iax12,coord2iax12             ! 2D target along axes 1,2
  real(wp), dimension(:), pointer :: coord1iax13,coord3iax13             ! 2D target along axes 1,3
  integer :: lc1i,lc2i,lc3i                                              ! dataset length along the 3 coordinate axes

  !! these are the input data arrays interpolated in time to the present (presuming we've called update/timeinterp
  real(wp), dimension(:), pointer :: data0Dinow
  real(wp), dimension(:,:), pointer :: data1Dax1inow
  real(wp), dimension(:,:), pointer :: data1Dax2inow,data1Dax3inow
  real(wp), dimension(:,:,:), pointer :: data2Dax23inow
  real(wp), dimension(:,:,:), pointer :: data2Dax12inow,data2dax13inow
  real(wp), dimension(:,:,:,:), pointer :: data3Dinow

  real(wp), dimension(2) :: tref                                     ! times for two input frames bracketting current time
  real(wp) :: tnow                                                   ! time corresponding to data in *now arrays, viz current time insofar as this object knows
  real(wp) :: dt                                                     ! time step for *this input data object*
  integer, dimension(3,2) :: ymdref                                     ! last dim is for prev,next
  real(wp), dimension(2) :: UTsecref

  contains
    !! top-level user-intended procedures; goal is to have only these called to interact with object
    procedure(initproc), deferred :: init        ! set up object for first time step:  call read in grid, set sizes, init_storage,
                                                 !   call prime_data, set data cadence based on some input
    procedure :: update                          ! check to see if new file needs to be read and read accordingly (will need to call deferred loaddata)
    procedure :: get_locationsi                  ! (no-op, extensions need to override) return a pointer to some locations to be used directly by user
    procedure :: get_datainow_ptr                ! (no-op) extensions need to return a pointer to a place where data can directly be fed
    procedure :: set_datainow                    ! (no-op, extensions shoudl override) user wants to directly set data from locations returned by get_locationsi

    !! internal/fine-grained control
    procedure :: set_sizes             ! initiate sizes for coordinate axes and number of datasets of different dimensionality
    procedure :: set_coords            ! fill interpolant coordinate arrays
    procedure :: set_cadence           ! fill dt variable for this instance of input
    procedure :: set_name              ! assign a character string name to our dataset
    procedure :: set_source            ! set the source directory for the input data
    procedure :: init_storage          ! wrapper routine to set up arrays once sizes are known/set
    procedure :: spaceinterp           ! interpolate spatially
    procedure :: timeinterp            ! interpolate in time based on data presently loaded into spatial arrays
    procedure :: dissociate_pointers   ! clear out memory and reset and allocation status flags
    procedure :: prime_data            ! load data buffers so that the object is ready for the first time step
    procedure :: update_simple         ! basic update steps needed for almost any datasets kind
    procedure :: init_storage_simple   ! basic allocations needed for most datasets

    !! internal, data kind specific
    procedure(coordisetproc), deferred :: set_coordsi       ! use grid data to compute coordinates of the interpolation sites
    procedure(loadgridsizeproc), deferred :: load_size          ! get size information from source data directory
    procedure(loadproc), deferred :: load_data              ! read data from file (possibly one array at a time) and spatially interpolate and store it in the appropriate arrays
    procedure(loadgridsizeproc), deferred :: load_grid              ! get grid information from source data directory
end type inputdata


!> interfaces for deferred procedures
abstract interface
  subroutine initproc(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    import inputdata,wp,curvmesh,gemini_cfg
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel
    real(wp), intent(in) :: dtdata
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
  end subroutine initproc

  subroutine coordisetproc(self,cfg,x)
    import inputdata,gemini_cfg,curvmesh
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in)  :: x
  end subroutine coordisetproc

  subroutine loadproc(self,t,dtmodel,ymdtmp,UTsectmp)
    import inputdata,wp
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
  end subroutine loadproc

  subroutine loadgridsizeproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine loadgridsizeproc
end interface

contains
  !> Load/store size variables.  By default this is going to get the sizes for the interpolated data from the grid.  If this
  !    is not the desired behavior then the type extension should override this procedure.
  subroutine set_sizes(self, &
                     l0D, &
                     l1Dax1,l1Dax2,l1Dax3, &
                     l2Dax23,l2Dax12,l2Dax13, &
                     l3D, &
                     x)
    class(inputdata), intent(inout) :: self
    integer, intent(in) :: l0D
    integer, intent(in) :: l1Dax1,l1Dax2,l1Dax3
    integer, intent(in) :: l2Dax23,l2Dax12,l2Dax13
    integer, intent(in) :: l3D
    class(curvmesh), intent(in) :: x        ! sizes for interpolation sites taken from grid

    ! Note:  these should be set by load_data deferred procedure
    ! coordinate axis sizes for input data
    !self%lc1=lc1; self%lc2=lc2; self%lc3=lc3;
    if (.not. self%flagdatasize) error stop 'inpudata:set_sizes() - must set input datasize first using load_size()'

    ! number of different types of data
    self%l0D=l0D
    self%l1Dax1=l1Dax1; self%l1Dax2=l1Dax2; self%l1Dax3=l1Dax3;
    self%l2Dax23=l2Dax23; self%l2Dax12=l2Dax12; self%l2Dax13=l2Dax13;
    self%l3D=l3D

    ! coordinate axis sizes for interpolation sites
    self%lc1i=x%lx1; self%lc2i=x%lx2; self%lc3i=x%lx3;

    ! check that the user is trying something sensible
    if (self%lc1==1 .and. self%lc1i/=1 .or. self%lc2==1 .and. self%lc2i/=1 &
            .or. self%lc3==1 .and. self%lc3i/=1) then
      if (self%flagforcenative) then
        print*, '  Warning:  native array rank forced for interpolations...'
      else if (self%flagallow2D3D) then
        print*, '  Warning:  allowing 2D to 3D spatial interpolations...'
      else
        print*, '  Dataset:  ',self%dataname,'  ',self%lc1,self%lc1i,self%lc2,self%lc2i,self%lc3,self%lc3i
        error stop 'inputdata:set_sizes() - singleton dimensions must be same for source and destination.'
      end if
    end if

    ! check dipole mesh or not
    select type (x)
      class is (dipolemesh)
        self%flagdipmesh=.true.
      class default
        self%flagdipmesh=.false.
    end select

    ! flag sizes as assigned
    self%flagsizes=.true.
  end subroutine set_sizes


  !> wrapper for allocations, should be overridden by type extensions if needed
  subroutine init_storage(self)
    class(inputdata), intent(inout) :: self

    call self%init_storage_simple()
  end subroutine init_storage


  !> allocate space to store inputdata
  subroutine init_storage_simple(self)
    class(inputdata), intent(inout) :: self
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

    ! NOTE: type extensions are reponsible for zeroing out any arrays they will use...

    ! input data coordinate arrays (presume plaid)
    allocate(self%coord1(lc1),self%coord2(lc2),self%coord3(lc3))

    ! interpolation site arrays (note these are flat, i.e. rank 1), if one needed to save space by not allocating unused block
    !   could override this procedure...
    allocate(self%coord1i(lc1i*lc2i*lc3i),self%coord2i(lc1i*lc2i*lc3i),self%coord3i(lc1i*lc2i*lc3i))
    ! coordinate sites for singleton axes depend on mangling of data
    if (self%flagdipmesh) then    ! mangle 2,3 sizes
      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc3i),self%coord3iax3(lc2i))
    else
      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
    end if

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
    allocate(self%data3Di(lc1i,lc2i,lc3i,l3D,2))

    ! allocate object arrays at interpolation sites for current time.  FIXME: do we even need to store permanently?
    allocate(self%data0Dinow(l0D))
    allocate(self%data1Dax1inow(lc1i,l1Dax1), self%data1Dax2inow(lc2i,l1Dax2), self%data1Dax3inow(lc3i,l1Dax3))
    allocate(self%data2Dax23inow(lc2i,lc3i,l2Dax23), self%data2Dax12inow(lc1i,lc2i,l2Dax12), self%data2Dax13inow(lc1i,lc3i,l2Dax13))
    allocate(self%data3Dinow(lc1i,lc2i,lc3i,l3D))

    self%flagalloc=.true.
  end subroutine init_storage_simple


  !> set the dataset names
  subroutine set_name(self,datasetstr)
    class(inputdata), intent(inout) :: self
    character(*), intent(in) :: datasetstr

    self%dataname=datasetstr
  end subroutine set_name


  !> set the cadence of the dataset
  subroutine set_cadence(self,dtdata)
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: dtdata

    self%dt=dtdata
    self%flagcadence=.true.
  end subroutine set_cadence


  !> set the source directory location
  subroutine set_source(self,sourcedir)
    class(inputdata), intent(inout) :: self
    character(*), intent(in) :: sourcedir

    self%sourcedir=sourcedir
    self%flagsource=.true.
  end subroutine set_source


  !> set the input data coordinates
  subroutine set_coords(self,c1,c2,c3)
    class(inputdata), intent(inout) :: self
    real(wp), dimension(:) :: c1,c2,c3

    if (.not. self%flagalloc) error stop ' inputdata:set_coords() - must allocate space prior to setting interpolant coordinates'

    self%coord1=c1
    self%coord2=c2
    self%coord3=c3
  end subroutine set_coords


  !> "prime" data at the beginning of the simulation so that proper inputs can be derived/interpolated for the first time step
  !     Note that we need to separate any activity that isn't directly related to input data (e.g. background states, etc.) from
  !     this routine so that it purely acts on properties of the inputdata class/type
  subroutine prime_data(self,cfg,x,dtmodel,ymd,UTsec)
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg      ! need start date, etc. of the simulation to prime
    class(curvmesh), intent(in) :: x         ! must have the grid to call update
    real(wp), intent(in) :: dtmodel             ! model time and time step
    integer, dimension(3), intent(in) :: ymd ! date from which we need to prime input data
    real(wp), intent(in) :: UTsec            ! time from which we need to prime input data
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp

    ! fIXME: unused variables

    !! do some really basic error checking
    if (self%flagdoinput) then
      if (.not. self%flagalloc) error stop 'inputdata:prime_data() - must allocate data arrays prior to priming'
      if (.not. self%flagcadence) error stop 'inputdata:prime_data() - must specify data cadence before priming'
      print*, '  Priming dataset:  ',self%dataname

      !! find the last input data preceding the milestone/initial condition that we start with
      !    The arguments here coorespond to start datetime of simulations, time of first step for this run (different
      !    if doing a restart) and then the tmp vars which are the time of the last input file.  dtdata is cadence.
      call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,self%dt,ymdtmp,UTsectmp)
      ! FIXME: just set to time of last neutral frame...
      !ymdtmp=ymd
      !UTsectmp=UTsec

      !! Loads the neutral input file corresponding to the "first" time step of the simulation to prevent the first interpolant
      !  from being zero and causing issues with restart simulations.  I.e. make sure the neutral buffers are primed for restart
      !  This requires us to load file input twice, once corresponding to the initial frame and once for the "first, next" frame.
      ! FIXME: need to keep self%ymd, etc. in sync?  Update will do this?  YES
      !self%tref(1)=UTsectmp-UTsec-2*self%dt
      !self%tref(2)=self%tref(1)+self%dt

      self%tref(1)=(UTsectmp-cfg%UTsec0)-2*self%dt
      self%tref(2)=self%tref(1)+self%dt

      !                         ' This is a workaround to insure compatibility with restarts...',ymdtmp,UTsectmp
      !! We essentially are loading up the data corresponding to halfway betwween -dtneu and t0 (zero).  This will load
      !   two time levels back so when tprev is incremented twice it will be the true tprev corresponding to first time step
      call self%update(cfg,dtmodel,self%tref(2)+self%dt/2,x,ymdtmp,UTsectmp-self%dt)  !abs time arg to be < 0

      !! Now compute perturbations for the present time (zero), this moves the primed variables in next into prev and then
      !  loads up a current state so that we get a proper interpolation for the first time step.
      !call self%update(cfg,dtmodel,0._wp,x,ymdtmp,UTsectmp)    !t-dt so we land exactly on start time
      call self%update(cfg,dtmodel,self%tref(2)+3/2*self%dt,x,ymdtmp,UTsectmp)    !t-dt so we land exactly on start time

      self%flagprimed=.true.

      print*, 'prime times:  ',self%tref(1),self%tref(2)
      print*, 'prime reference date:  ', self%ymdref(:,1),self%UTsecref(1),self%ymdref(:,2),self%UTsecref(2)
    end if
  end subroutine prime_data


  !> default wrapper for updates, extension can override as needed or just adopt this default
  subroutine update(self,cfg,dtmodel,t,x,ymd,UTsec)
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel    ! need both model and input data time stepping
    real(wp), intent(in) :: t                 ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x         ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations

    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)
  end subroutine update


  !> routine to execute various steps needed to update input data to present time
  subroutine update_simple(self,cfg,dtmodel,t,x,ymd,UTsec)
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel    ! need both model and input data time stepping
    real(wp), intent(in) :: t                 ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x         ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we compute perturbations

    integer, dimension(3) :: ymdtmp          ! these hold the incremented date following reading of new file
    real(wp) :: UTsectmp

    !! basic error checking
    if (.not. self%flagalloc) error stop 'inputdata:update() - must allocate array space prior to update...'
    if (.not. self%flagcadence) error stop 'inputdata:update() - must define cadence first...'
    if (.not. self%flagsource) error stop 'inputdata:update() - must define source data directory...'

    !print*, 'entering update',ymd,UTsec,t,t+dtmodel/2,self%tref(1),self%tref(2)
    !print*, '    ',self%ymdref(:,1),self%UTsecref(1),self%ymdref(:,2),self%UTsecref(2)

    !! see if we need to load new data into the buffer; negative time means that we need to load the first frame
    if (t+dtmodel/2 >= self%tref(2) .or. t < 0) then
      !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
      if (self%flagfirst) then
        !initialize dates
        self%ymdref(:,1)=ymd
        self%UTsecref(1)=UTsec
        self%ymdref(:,2)=self%ymdref(:,1)
        self%UTsecref(2)=self%UTsecref(1)
        self%flagfirst=.false.
      end if
      if (.not. self%flagcoordsi) then     !means this is the first time we've tried to load neutral simulation data, should we check for a previous neutral file to load??? or just assume everything starts at zero?  This needs to somehow check for an existing file under certain conditiosn, maybe if it==1???  Actually we don't even need that we can just check that the neutral grid is allocated (or not)
        !Create a neutral grid, do some allocations and projections
        call self%set_coordsi(cfg,x)    ! cfg needed to convey optional parameters about how the input data are to be
                                       !interpreted
      end if

      !Read in neutral data from a file
      call self%load_data(t,dtmodel,ymdtmp,UTsectmp)

      !Spatial interpolation for the frame we just read in
      call self%spaceinterp()

      !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
      self%tref(1)=self%tref(2)
      self%UTsecref(1)=self%UTsecref(2)
      self%ymdref(:,1)=self%ymdref(:,2)

      self%tref(2)=self%tref(1)+self%dt
      self%UTsecref(2)=UTsectmp
      self%ymdref(:,2)=ymdtmp
    end if !done loading frame data...

    !Interpolation in time
    call self%timeinterp(t,dtmodel)
  end subroutine update_simple


  !> use data stored in input arrays to interpolate onto grid sites for "next" dataset.  There may be a need here to
  !    accommodate singleton dimension naturally to void having to define extensions for different types of interp...
  subroutine spaceinterp(self)
    class(inputdata),intent(inout) :: self
    integer :: iparm
    real(wp), dimension(:), allocatable :: tempdata
    integer :: lc1i,lc2i,lc3i,lc1,lc2,lc3
    real(wp), dimension(:), pointer :: coord1,coord2,coord3
    real(wp), dimension(:), pointer :: coord1i,coord2i,coord3i
    real(wp), dimension(:), pointer :: coord1iax1,coord2iax2,coord3iax3
    real(wp), dimension(:), pointer :: coord2iax23,coord3iax23
    real(wp), dimension(:), pointer :: coord1iax12,coord2iax12
    real(wp), dimension(:), pointer :: coord1iax13,coord3iax13

    ! FIXME: possibly needs some more error checking

    ! for convenience
    lc1i=self%lc1i; lc2i=self%lc2i; lc3i=self%lc3i;
    lc1=self%lc1; lc2=self%lc2; lc3=self%lc3;
    coord1=>self%coord1; coord2=>self%coord2; coord3=>self%coord3;
    coord1i=>self%coord1i; coord2i=>self%coord2i; coord3i=>self%coord3i;
    coord1iax1=>self%coord1iax1; coord2iax2=>self%coord2iax2; coord3iax3=>self%coord3iax3;
    coord2iax23=>self%coord2iax23; coord3iax23=>self%coord3iax23;
    coord1iax12=>self%coord1iax12; coord2iax12=>self%coord2iax12;
    coord1iax13=>self%coord1iax13; coord3iax13=>self%coord3iax13;

    !> 1D arrays varying along the 1-axis
    if (self%l1Dax1>0) then
      self%data1Dax1i(:,:,1)=self%data1Dax1i(:,:,2)     ! save old data!!!
      allocate(tempdata(self%lc1i))
      do iparm=1,self%l1Dax1
        tempdata(:)=interp1(coord1,self%data1Dax1(:,iparm),coord1iax1)
        self%data1Dax1i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 2-axis
    if (self%l1Dax2>0) then
      self%data1Dax2i(:,:,1)=self%data1Dax2i(:,:,2)
      allocate(tempdata(self%lc2i))
      do iparm=1,self%l1Dax2
        tempdata(:)=interp1(coord2,self%data1Dax2(:,iparm),coord2iax2)
        self%data1Dax2i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 3-axis
    if (self%l1Dax3>0) then
      self%data1Dax3i(:,:,1)=self%data1Dax3i(:,:,2)
      allocate(tempdata(self%lc3i))
      do iparm=1,self%l1Dax3
        tempdata(:)=interp1(coord3,self%data1Dax3(:,iparm),coord3iax3)
        self%data1Dax3i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 2,3 axes; be sure to check singleton axes and change interp shape accordingly
    if (self%l2Dax23>0) then
      self%data2Dax23i(:,:,:,1)=self%data2Dax23i(:,:,:,2)
      if (lc2>1 .and. lc3>1 .or. self%flagforcenative) then    ! normal 2D dataset
        allocate(tempdata(self%lc2i*self%lc3i))
        do iparm=1,self%l2Dax23
          tempdata(:)=interp2(coord2,coord3,self%data2Dax23(:,:,iparm),coord2iax23,coord3iax23)
          self%data2Dax23i(:,:,iparm,2)=reshape(tempdata,[lc2i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc2>1 .and. lc3==1) then
        allocate(tempdata(self%lc2i*self%lc3i))   ! in case lc2,3i are mangled just allociate using total size...
        do iparm=1,self%l2Dax23
          tempdata(:)=interp1(coord2,self%data2Dax23(:,1,iparm),coord2iax23)
          self%data2Dax23i(:,:,iparm,2)=reshape(tempdata,[lc2i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc2==1 .and. lc3>1) then
        allocate(tempdata(self%lc2i*self%lc3i))
        do iparm=1,self%l2Dax23
          tempdata(:)=interp1(coord3,self%data2Dax23(1,:,iparm),coord3iax23)
          self%data2Dax23i(:,:,iparm,2)=reshape(tempdata,[lc2i,lc3i])
        end do
        deallocate(tempdata)
      else
        error stop 'inputdata:spaceinterp() - cannot determine type of interpolation for data2Dax23'
      end if
    end if

    !> 2D arrays varying along the 1,2 axes
    if (self%l2Dax12>0) then
      self%data2Dax12i(:,:,:,1)=self%data2Dax12i(:,:,:,2)
      if (lc1>1 .and. lc2>1 .or. self%flagforcenative) then
        allocate(tempdata(self%lc1i*self%lc2i))
        do iparm=1,self%l2Dax12
          tempdata(:)=interp2(coord1,coord2,self%data2Dax12(:,:,iparm),coord1iax12,coord2iax12)
          self%data2Dax12i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i])
        end do
        deallocate(tempdata)
      else if (lc1>1 .and. lc2==1) then
        allocate(tempdata(self%lc1i*self%lc2i))
        do iparm=1,self%l2Dax12
          tempdata(:)=interp1(coord1,self%data2Dax12(:,1,iparm),coord1iax12)
          self%data2Dax12i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i])
        end do
        deallocate(tempdata)
      else if (lc1==1 .and. lc2>1) then
        allocate(tempdata(self%lc1i*self%lc2i))
        do iparm=1,self%l2Dax12
          tempdata(:)=interp1(coord2,self%data2Dax12(1,:,iparm),coord2iax12)
          self%data2Dax12i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i])
        end do
        deallocate(tempdata)
      else
        error stop 'inputdata:spaceinterp() - cannot determine type of interpolation for data2Dax12'
      end if
    end if

    !> 2D arrays varying along the 1,3 axes
    if (self%l2Dax13>0) then
      self%data2Dax13i(:,:,:,1)=self%data2Dax13i(:,:,:,2)
      if (lc1>1 .and. lc3>1 .or. self%flagforcenative) then
        allocate(tempdata(self%lc1i*self%lc3i))
        do iparm=1,self%l2Dax13
          tempdata(:)=interp2(coord1,coord3,self%data2Dax13(:,:,iparm),coord1iax13,coord3iax13)
          self%data2Dax13i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc1>1 .and. lc3==1) then
        allocate(tempdata(self%lc1i*self%lc3i))
        do iparm=1,self%l2Dax13
          tempdata(:)=interp1(coord1,self%data2Dax13(:,1,iparm),coord1iax13)
          self%data2Dax13i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc1==1 .and. lc3>1) then
        allocate(tempdata(self%lc1i*self%lc3i))
        do iparm=1,self%l2Dax13
          tempdata(:)=interp1(coord3,self%data2Dax13(1,:,iparm),coord3iax13)
          self%data2Dax13i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc3i])
        end do
        deallocate(tempdata)
      else
        error stop 'inputdata:spaceinterp() - cannot determine type of interpolation for data2Dax13'
      end if
    end if

    !> 3D arrays varying along all axes, check for singleton axes...
    if (self%l3D>0) then
      self%data3Di(:,:,:,:,1)=self%data3Di(:,:,:,:,2)
      if (lc1>1 .and. lc2>1 .and. lc3>1 .or. self%flagforcenative) then    ! forcenative because sometimes we interp 2D->3D, e.g. for neutral axisymmetric inputs
        allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
        do iparm=1,self%l3D
          tempdata(:)=interp3(coord1,coord2,coord3,self%data3D(:,:,:,iparm),coord1i,coord2i,coord3i)
          self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc1>1 .and. lc2>1 .and. lc3==1) then
        allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
        do iparm=1,self%l3D
          tempdata(:)=interp2(coord1,coord2,self%data3D(:,:,1,iparm),coord1i,coord2i)
          self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc1>1 .and. lc2==1 .and. lc3>1) then
        allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
        do iparm=1,self%l3D
          tempdata(:)=interp2(coord1,coord3,self%data3D(:,1,:,iparm),coord1i,coord3i)
          self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
        end do
        deallocate(tempdata)
      else if (lc1==1 .and. lc2>1 .and. lc3>1) then
        allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
        do iparm=1,self%l3D
          tempdata(:)=interp2(coord2,coord3,self%data3D(1,:,:,iparm),coord2i,coord3i)
          self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
        end do
        deallocate(tempdata)
      else
        error stop 'inputdata:spaceinterp() - cannot determine type of interpolation for data3D'
      end if
    end if
  end subroutine spaceinterp


  !> interpolate data in time (requires first loading and interpolating in space)
  subroutine timeinterp(self,t,dt)
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: t,dt
    real(wp) :: slope
    integer :: ic1,ic2,ic3,iparm
    integer :: lc1i,lc2i,lc3i

    ! convenience vars
    lc1i=self%lc1i; lc2i=self%lc2i; lc3i=self%lc3i;

    !> interpolate scalars in time, never enter loop for params of 0 size
    do iparm=1,self%l0D
      slope=(self%data0Di(iparm,2)-self%data0Di(iparm,1))/(self%tref(2)-self%tref(1))
      self%data0Dinow(iparm)=self%data0Di(iparm,1)+slope*(t+dt/2 -self%tref(1))
    end do

    !> 1D arrays varying along the 1-axis
    do iparm=1,self%l1Dax1
      do ic1=1,lc1i
        slope=(self%data1Dax1i(ic1,iparm,2)-self%data1Dax1i(ic1,iparm,1))/(self%tref(2)-self%tref(1))
        self%data1Dax1inow(ic1,iparm)=self%data1Dax1i(ic1,iparm,1)+slope*(t+dt/2 -self%tref(1))
      end do
    end do

    !> 1D arrays varying along the 2-axis
    do iparm=1,self%l1Dax2
      do ic2=1,lc2i
        slope=(self%data1Dax2i(ic2,iparm,2)-self%data1Dax2i(ic2,iparm,1))/(self%tref(2)-self%tref(1))
        self%data1Dax2inow(ic2,iparm)=self%data1Dax2i(ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
      end do
    end do

    !> 1D arrays varying along the 3-axis
    do iparm=1,self%l1Dax3
      do ic3=1,lc3i
        slope=(self%data1Dax3i(ic3,iparm,2)-self%data1Dax3i(ic3,iparm,1))/(self%tref(2)-self%tref(1))
        self%data1Dax3inow(ic3,iparm)=self%data1Dax3i(ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
      end do
    end do

    !> 2D arrays varying along the 2,3 axes
    do iparm=1,self%l2Dax23
      do ic3=1,lc3i
        do ic2=1,lc2i
          slope=(self%data2Dax23i(ic2,ic3,iparm,2)-self%data2Dax23i(ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
          self%data2Dax23inow(ic2,ic3,iparm)=self%data2Dax23i(ic2,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end do

    !> 2D arrays varying along the 1,2 axes
    do iparm=1,self%l2Dax12
      do ic2=1,lc2i
        do ic1=1,lc1i
          slope=(self%data2Dax12i(ic1,ic2,iparm,2)-self%data2Dax12i(ic1,ic2,iparm,1))/(self%tref(2)-self%tref(1))
          self%data2Dax12inow(ic1,ic2,iparm)=self%data2Dax12i(ic1,ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end do

    !> 2D arrays varying along the 1,3 axes
    do iparm=1,self%l2Dax13
      do ic3=1,lc3i
        do ic1=1,lc1i
          slope=(self%data2Dax13i(ic1,ic3,iparm,2)-self%data2Dax13i(ic1,ic3,iparm,1))/(self%tref(2)-self%tref(1))
          self%data2Dax13inow(ic1,ic3,iparm)=self%data2Dax13i(ic1,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end do

    !> 3D arrays varying along all axes (obv.)
    do iparm=1,self%l3D
      do ic3=1,lc3i
        do ic2=1,lc2i
          do ic1=1,lc1i
            slope=(self%data3Di(ic1,ic2,ic3,iparm,2)-self%data3Di(ic1,ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
            self%data3Dinow(ic1,ic2,ic3,iparm)=self%data3Di(ic1,ic2,ic3,iparm,1)+slope*(t+dt/2-self%tref(1))
          end do
        end do
      end do
    end do

    !> update current time in object
    self%tnow=t+dt/2
  end subroutine timeinterp


  !> These will do nothing for now, can override with custom code as needed
  subroutine get_locationsi(self,flagallpts,zlims,xlims,ylims,zvals,xvals,yvals,datavals)
    class(inputdata), intent(inout) :: self
    logical, intent(in) :: flagallpts
    real(wp), dimension(2), intent(in) :: zlims,xlims,ylims    ! global boundary of neutral grid we are accepting data from
    real(wp), dimension(:), pointer, intent(inout) :: zvals,xvals,yvals
    real(wp), dimension(:,:), pointer, intent(inout) :: datavals

    print*, 'WARNING:  triggered no-op get_locationsi, use an extension with a full implementation'
    return
  end subroutine


  function get_datainow_ptr(self) result(datavals)
    class(inputdata), intent(inout) :: self
    real(wp), dimension(:,:), pointer :: datavals 

    print*, 'WARNING:  triggered no-op get_locationsi, use an extension with a full implementation'
    return
  end function get_datainow_ptr


  !> We assume that the get_locationsi will provide a memory space for the results which are stored in the object extension
  !    so no additional inputs are needed to copy those data out into the proper object arrays.  
  subroutine set_datainow(self)
    class(inputdata), intent(inout) :: self

    print*, 'WARNING:  triggered no-op set_datainow, use an extension with a full implementation'           
    return
  end subroutine set_datainow


  !> deallocate memory and dissociated pointers for generic array data
  subroutine dissociate_pointers(self)
    class(inputdata), intent(inout) :: self

    if (self%flagalloc) then
      deallocate(self%data0D)
      deallocate(self%data1Dax1, self%data1Dax2, self%data1Dax3)
      deallocate(self%data2Dax23, self%data2Dax12, self%data2Dax13)
      deallocate(self%data3D)

      deallocate(self%data0Di)
      deallocate(self%data1Dax1i, self%data1Dax2i, self%data1Dax3i)
      deallocate(self%data2Dax23i, self%data2Dax12i, self%data2Dax13i)
      deallocate(self%data3Di)

      deallocate(self%coord1,self%coord2,self%coord3)
      deallocate(self%coord1i,self%coord2i,self%coord3i)
    end if

    self%flagalloc=.false.
    self%flagprimed=.false.
    self%flagcoordsi=.false.
  end subroutine dissociate_pointers
end module inputdataobj
