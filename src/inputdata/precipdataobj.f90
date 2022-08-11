module precipdataobj

! type extension for file-based precipitation data input.  Assumes parallel communication between root/workers for data
! distribution.

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,debug,pi
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize2,get_grid2,get_precip
use timeutils, only: dateinc,date_filename

implicit none (type, external)
private
public :: precipdata

type, extends(inputdata) :: precipdata
  ! coordinate for input precipitation data, and storage
  real(wp), dimension(:), pointer :: mlonp,mlatp
  integer, pointer :: llon,llat
  real(wp), dimension(:,:), pointer :: Qp,E0p
  real(wp), dimension(:,:), pointer :: Qpiprev,E0piprev
  real(wp), dimension(:,:), pointer :: Qpinext,E0pinext
  real(wp), dimension(:,:), pointer :: Qpinow,E0pinow

  contains
    ! overriding procedures
    procedure :: set_sizes=>set_sizes_precip

    ! deferred bindings
    procedure :: init=>init_precip
    procedure :: set_coordsi=>set_coordsi_precip
    procedure :: load_data=>load_data_precip
    procedure :: load_grid=>load_grid_precip
    procedure :: load_size=>load_size_precip     ! load the size of the input data files
    final :: destructor
end type precipdata

contains
  !> need to override set_sizes so to account for fact that target interpolation is to a 2D grid, by default object will
  !    assume 3D and get the sizes from the simulation grid
  subroutine set_sizes_precip(self, &
                     l0D, &
                     l1Dax1,l1Dax2,l1Dax3, &
                     l2Dax23,l2Dax12,l2Dax13, &
                     l3D, &
                     x)
    class(precipdata), intent(inout) :: self
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
    ! coordinate axis sizes for interpolation states
    !select type (x)
    !  class is (dipolemesh)
    !    print*, ' precipdata:  detected dipole mesh...'
    !    self%lc1i=x%lx1;       ! note this dataset has 1D and 2D target interpolation grid
    !    self%lc2i=x%lx3; self%lc3i=x%lx2;    ! dipolemesh mesh permuted ~alt,lat,lon  more or less...
    !  class default
        self%lc1i=x%lx1;       ! note this dataset has 1D and 2D target interpolation grid
        self%lc2i=x%lx2; self%lc3i=x%lx3;
    !end select

    ! check that the user is trying something sensible
    !if (self%lc1==1 .and. self%lc1i/=1 .or. self%lc2==1 .and. self%lc2i/=1 &
    !        .or. self%lc3==1 .and. self%lc3i/=1) then
    !  error stop 'inputdata:set_sizes() - singleton dimensions must be same for source and destination.'
    !end if

    ! flag sizes as assigned
    self%flagsizes=.true.
  end subroutine set_sizes_precip


  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_precip(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(precipdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character(*), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtmodel,dtdata                      ! model time step and cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='electron precipitation'
    call self%set_name(strname)
    self%flagdoinput=cfg%flagprecfile/=0

    ! read the simulation size from the source directory and allocate arrays
    allocate(self%lc1,self%lc2,self%lc3)      ! these are pointers
    self%llon=>self%lc2; self%llat=>self%lc3;
    call self%load_size()
    call self%set_sizes(0, &
                       0,0,0, &
                       2,0,0, &
                       0, &
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%mlonp=>self%coord2; self%mlatp=>self%coord3;
    call self%load_grid()

    ! set input data array pointers to faciliate easy to read input code; these may or may not be helpful to user
    self%Qp=>self%data2Dax23(:,:,1)
    self%E0p=>self%data2Dax23(:,:,2)
    self%Qpiprev=>self%data2Dax23i(:,:,1,1)
    self%Qpinext=>self%data2Dax23i(:,:,1,2)
    self%E0piprev=>self%data2Dax23i(:,:,2,1)
    self%E0pinext=>self%data2Dax23i(:,:,2,2)
    self%Qpinow=>self%data2Dax23inow(:,:,1)
    self%E0pinow=>self%data2Dax23inow(:,:,2)

    ! must initialize prev state or else the first set of data will not be interpolated correctly
    self%Qpiprev=0.0
    self%Qpinext=0.0
    self%E0piprev=100.0
    self%E0pinext=100.0

    ! set to start time of simulation - will be set first time update is called
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_precip


  !> get the input grid size from file, all workers will just call this sicne this is a one-time thing
  subroutine load_size_precip(self)
    class(precipdata), intent(inout) :: self

    ! basic error checking
    if (.not. self%flagsource) error stop 'precipdata:load_size_precip() - must define a source directory first'

    ! read sizes
    print '(/,A,/,A)', 'Precipitation input:','--------------------'
    print '(A)', 'READ precipitation size from: ' // self%sourcedir
    call get_simsize2(self%sourcedir // "/simsize.h5", llon=self%llon, llat=self%llat)

    print '(A,2I6)', 'Precipitation size: llon,llat:  ',self%llon,self%llat
    if (self%llon < 1 .or. self%llat < 1) then
     print*, '  precipitation grid size must be strictly positive: ' //  self%sourcedir
     error stop
    end if

    ! set dim 1 size to null since inherently not used
    self%lc1=0

    ! flag to denote input data size is set
    self%flagdatasize=.true.
  end subroutine load_size_precip


  !> get the grid information from a file, all workers will just call this since one-time
  subroutine load_grid_precip(self)
    class(precipdata), intent(inout) :: self

    ! read grid data
    call get_grid2(self%sourcedir // "/simgrid.h5", self%mlonp, self%mlatp)

    print '(A,4F9.3)', 'Precipitation mlon,mlat extent:  ',minval(self%mlonp(:)),maxval(self%mlonp(:)), &
                                                           minval(self%mlatp(:)),maxval(self%mlatp(:))
    if(.not. all(ieee_is_finite(self%mlonp))) error stop 'precipBCs_fileinput: mlon must be finite'
    if(.not. all(ieee_is_finite(self%mlatp))) error stop 'precipBCs_fileinput: mlat must be finite'
  end subroutine load_grid_precip


  !> set target coordinates for interpolation sights
  subroutine set_coordsi_precip(self,cfg,x)
    class(precipdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg     ! presently not used but possibly eventually?
    class(curvmesh), intent(in) :: x
    integer :: ix2,ix3,iflat

    iflat = cfg%potsolve
    !! avoid unused argument warning

    ! set full 2D target coordinates along axes 2,3 - these are the only targets we have for precipitation data
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        iflat=(ix3-1)*x%lx2+ix2
        self%coord2iax23(iflat)=x%phi(x%lx1,ix2,ix3)*180._wp/pi
        self%coord3iax23(iflat)=90._wp - x%theta(x%lx1,ix2,ix3)*180._wp/pi
      end do
    end do

    self%flagcoordsi=.true.
  end subroutine set_coordsi_precip


  !> have all processes read in data from file to avoid any message passing
  subroutine load_data_precip(self,t,dtmodel,ymdtmp,UTsectmp)
    class(precipdata), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warnings

    !! all workers should update the date
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt, ymdtmp, UTsectmp)

    !! this read must be done repeatedly through simulation so have only root do file io
    !print*, '  date and time:  ',ymdtmp,UTsectmp
    !print*, '  precip filename:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp)
    ! read in the data for the "next" frame from file
    call get_precip(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", self%Qp, self%E0p)

    !print*, ' precip data succesfully input...'
  end subroutine load_data_precip


  !> destructor needs to clear memory out
  subroutine destructor(self)
    type(precipdata), intent(inout) :: self

    call self%dissociate_pointers()
  end subroutine destructor
end module precipdataobj
