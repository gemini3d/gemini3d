module precipdataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata,initproc,coordsisetproc,loadproc
use meshobj, only: curvmesh
use config, only: gemini_cfg
use reader, only: get_simsize2,get_grid2,get_precip

implicit none (type, external)

type, extends(inputdata) :: precipdata
  ! coordinate for input precipitation data, and storage
  real(wp), dimension(:), pointer :: mlonp,mlatp
  integer, pointer :: llon,llat
  real(wp), dimension(:,:), pointer :: Qp,E0p

  contains
    procedure(initproc) :: init=>init_precip
    procedure(coordsisetproc) :: set_coordsi=>set_coordsi_precip
    procedure(loadproc) :: load_data=>load_data_precip
    procedure(loadproc) :: load_grid=>load_grid_precip
    procedure(loadproc) :: load_size=>load_size_precip
    final :: destructor
end type precipdata

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_precip(self,cfg,sourcedir,x,dtdata,ymd,UTsec)
    class(precipdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character, dimension(:), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtdata                      ! cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation

    ! tell our object where its data are and give the dataset a name
    self%set_source(sourcedir)
    self%set_name('electron precipitation')

    ! read the simulation size from the source directory and allocate arrays
    self%llon=>lc2; self%llat=>lc3
    call self%load_size()
    call self%set_sizes(lc1,lc2,lc3, &
                       0, &
                       0,0,0 &
                       2,0,0 &
                       0, &
                       x )
    call self%init_storage()
    self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%mlonp=>self%coord2; self%mlatp=>self%coord3;
    call self%load_grid()

    ! set input data array pointers to faciliate easy to read input code.  It's important to note that
    !   
    if (llon>1 .and. llat>1) then      ! full 2D input arrays
      Qp=>data2Dax23(:,:,1)
      E0p=>data2Dax23(:,:,2)
      Qpiprev=>data2Dax23i(:,:,1,1)
      Qpinext=>data2Dax23i(:,:,1,2)
      E0iprev=>data2Dax23i(:,:,2,1)
      E0inext=>data2Dax23i(:,:,2,2)
      Qpinow=>data2Dax23inow(:,:,1)
      E0inow=>data2Dax23inow(:,:,2)
    else if (llon==1 .and. llat>1) then    ! alt/lat simulation, data vary along 3rd axis
      Qp=>data1Dax3(:,1)
      E0p=>data1Dax3(:,2)
      Qpiprev=>data1Dax3i(:,1,1)
      Qpinext=>data1Dax3i(:,1,2)
      E0iprev=>data1Dax3i(:,2,1)
      E0inext=>data1Dax3i(:,2,2)
      Qpinow=>data1Dax3inow(:,1)
      E0inow=>data1Dax3inow(:,2)
    else if (llat==1 .and. llon>1) then    ! alt/lon simulation, input data vary along 2nd axis
      Qp=>data1Dax2(:,1)
      E0p=>data1Dax2(:,2)
      Qpiprev=>data1Dax2i(:,1,1)
      Qpinext=>data1Dax2i(:,1,2)
      E0iprev=>data1Dax2i(:,2,1)
      E0inext=>data1Dax2i(:,2,2)
      Qpinow=>data1Dax2inow(:,1)
      E0inow=>data1Dax2inow(:,2)
    else
      error stop 'precipdata:init_precip() - unsupported input file size configuration'
    end if

    ! prime input data
    call self%prime_data(cfg,x,ymd,UTsec)
  end subroutine init_precip


  !> get the input grid size from file, all workers will just call this
  subroutine load_size_precip(self)
    class(precipdata), intent(inout) :: self

    ! basic error checking
    if (.not. flagsource) error stop 'precipdata:load_size_precip() - must define a source directory first'

    ! read sizes
    print '(/,A,/,A)', 'Precipitation input:','--------------------'
    print '(A)', 'READ precipitation size from: ' // self%sourcedir
    call get_simsize2(sourcedir, llon=self%llon, llat=self%llat)

    print '(A,2I6)', 'Precipitation size: llon,llat:  ',llon,llat
    if (self%llon < 1 .or. self%llat < 1) error stop 'precipitation grid size must be strictly positive: ' //  self%sourcedir
  end subroutine load_size_precip


  !> get the grid information from a file, all workers will just call this
  subroutine load_grid_precip(self)
    class(precipdata), intent(inout) :: self

    ! read grid data
    call get_grid2(self%sourcedir, self%mlonp, self%mlatp)

    print '(A,4F9.3)', 'Precipitation mlon,mlat extent:  ',minval(self%mlonp(:)),maxval(self%mlonp(:)), &
                                                           minval(self%mlatp(:)),maxval(self%mlatp(:))
    if(.not. all(ieee_is_finite(self%mlonp))) error stop 'precipBCs_fileinput: mlon must be finite'
    if(.not. all(ieee_is_finite(self%mlatp))) error stop 'precipBCs_fileinput: mlat must be finite'
  end subroutine load_grid_precip

  subroutine load_data_precip(self)
    class(precipdata), intent(inout) :: self

    call get_precip(date_filename(self%sourcedir,ymd,UTsec))
  end subroutine load_data_precip
end module precipdataobj
