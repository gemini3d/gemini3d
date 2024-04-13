module solfluxdataobj

! type extension for file-based precipitation data input.  Assumes parallel communication between root/workers for data
! distribution.

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,debug,pi
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize2,get_grid2!,get_solflux
use timeutils, only: dateinc,date_filename

implicit none (type, external)
private
public :: solfluxdata

type, extends(inputdata) :: solfluxdata
  ! coordinates for input precipitation data, and storage
  real(wp), dimension(:), pointer :: glonp,glatp
  integer, pointer :: llon,llat,lalt
  real(wp), dimension(:,:,:), pointer :: Iinfp
  real(wp), dimension(:,:,:,:), pointer :: Iinfiprev,Iinfinext,Iinfinow

  ! work and target coordinates
  real(wp), dimension(:,:,:), allocatable :: glonimat,glatimat
  real(wp), dimension(:), pointer :: gloni,glati

  contains
    ! deferred bindings
    procedure :: init=>init_solflux
    procedure :: set_coordsi=>set_coordsi_solflux
    procedure :: load_data=>load_data_solflux
    procedure :: load_grid=>load_grid_solflux
    procedure :: load_size=>load_size_solflux     ! load the size of the input data files
    final :: destructor
end type solfluxdata

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_solflux(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(solfluxdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character(*), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtmodel,dtdata                      ! model time step and cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname

    ! need to allow interpolation from 2D to 3D
    self%flagallow2D3D=.true.

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='solar flux'
    call self%set_name(strname)
    !self%flagdoinput=cfg%flagsolfluxfile/=0

    ! read the simulation size from the source directory and allocate arrays
    allocate(self%lc1,self%lc2,self%lc3)      ! these are pointers
    self%llon=>self%lc2;  self%llat=>self%lc3;  self%lalt=>self%lc1;
    call self%load_size()
    call self%set_sizes(0, &
                       0,0,0, &
                       0,0,0, &     ! 22 different wavelength bins to interpolate for GEMINI's solar flux calculations
                       22, &
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%glonp=>self%coord2; self%glatp=>self%coord3;
    call self%load_grid()

    ! Set input data array pointers to faciliate easy to read input code; these may or may not be helpful to user
    !   We're treating solar flux as 2D for purposes of interpolation sources; however, convenient here to alias as a 3D array 
    !   with the 3rd axis being wavelength.  The target grid for interpolation is effectively glon,glat but treated
    !   as a 3D grid since we need a value of Iinf for every possible glon,glat on the grid and these are not plaid
    !   and cannot be mapped from a 2D array of glon,glat easily.  
    self%Iinfp=>self%data3D(1,:,:,:)
    self%Iinfiprev=>self%data3Di(:,:,:,:,1)
    self%Iinfinext=>self%data3Di(:,:,:,:,2)
    self%Iinfinow=>self%data3Dinow(:,:,:,:)

    ! must initialize prev state or else the first set of data will not be interpolated correctly
    self%Iinfiprev=0.0
    self%Iinfinext=0.0

    ! set to start time of simulation - will be set first time update is called
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_solflux


  !> get the input grid size from file, all workers will just call this sicne this is a one-time thing
  subroutine load_size_solflux(self)
    class(solfluxdata), intent(inout) :: self

    ! basic error checking
    if (.not. self%flagsource) error stop 'solfluxdata:load_size_solflux() - must define a source directory first'

    ! read sizes
    print '(/,A,/,A)', 'solflux input:','--------------------'
    print '(A)', 'READ solflux size from: ' // self%sourcedir
    call get_simsize2(self%sourcedir // "/simsize.h5", llon=self%llon, llat=self%llat)

    print '(A,2I6)', 'soflux size: llon,llat:  ',self%llon,self%llat
    if (self%llon < 1 .or. self%llat < 1) then
     print*, '  solflux grid size must be strictly positive: ' //  self%sourcedir
     error stop
    end if

    ! set dim 1 so we can use interpolation into a 3D array
    self%lc1=1

    ! flag to denote input data size is set
    self%flagdatasize=.true.
  end subroutine load_size_solflux


  !> get the grid information from a file, all workers will just call this since one-time
  subroutine load_grid_solflux(self)
    class(solfluxdata), intent(inout) :: self

    ! read grid data
    call get_grid2(self%sourcedir // "/simgrid.h5", self%glonp, self%glatp)

    print '(A,4F9.3)', 'Precipitation glon,glat extent:  ',minval(self%glonp(:)),maxval(self%glonp(:)), &
                                                           minval(self%glatp(:)),maxval(self%glatp(:))
    if(.not. all(ieee_is_finite(self%glonp))) error stop 'solfluxBCs_fileinput: glon must be finite'
    if(.not. all(ieee_is_finite(self%glatp))) error stop 'solfluxBCs_fileinput: glat must be finite'
  end subroutine load_grid_solflux


  !> set target coordinates for interpolation sights
  subroutine set_coordsi_solflux(self,cfg,x)
    class(solfluxdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg     ! presently not used but possibly eventually?
    class(curvmesh), intent(in) :: x
    integer :: ix1,ix2,ix3

    ! aliases for target interpolation sites
    self%gloni=>self%coord2i
    self%glati=>self%coord3i

    allocate(self%glonimat(1:x%lx1,1:x%lx2,1:x%lx3))     ! why not local variables?  FIXME
    allocate(self%glatimat,mold=self%glonimat)

    ! Target coordinates are 3D in this case...
    ! set full 2D target coordinates along axes 2,3 - these are the only targets we have for precipitation data
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          self%glonimat(ix1,ix2,ix3)=x%glon(ix1,ix2,ix3)
          self%glatimat(ix1,ix2,ix3)=x%glat(ix1,ix2,ix3)
        end do
      end do
    end do
    self%gloni=pack(self%glonimat,.true.)
    self%glati=pack(self%glatimat,.true.)

    deallocate(self%glonimat,self%glatimat)
    self%flagcoordsi=.true.
  end subroutine set_coordsi_solflux


  !> have all processes read in data from file to avoid any message passing
  subroutine load_data_solflux(self,t,dtmodel,ymdtmp,UTsectmp)
    class(solfluxdata), intent(inout) :: self
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
    !print*, '  precip filename:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp)
    ! read in the data for the "next" frame from file
    !call get_solflux(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", self%Qp, self%E0p)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !print*, ' precip data succesfully input...'
  end subroutine load_data_solflux


  !> destructor needs to clear memory out
  subroutine destructor(self)
    type(solfluxdata), intent(inout) :: self

    call self%dissociate_pointers()
  end subroutine destructor
end module solfluxdataobj
