module efielddataobj

! type extension for file-based electric field data input.  Assumes parallel communication between root/workers for data
! distribution.  An interesting aspect of this object is that it is/may never be used by anyone but root who initiates the
! calls to MUMPS and boundary conditions.

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize2,get_grid2,get_efield
!! note that only root uses this data object since this concerns the potential solver
use timeutils, only: dateinc,date_filename
use grid, only: lx1,lx2,lx2all,lx3,lx3all,gridflag

implicit none (type, external)
private
public :: efielddata

type, extends(inputdata) :: efielddata
  ! coordinate for input precipitation data, and storage
  real(wp), dimension(:), pointer :: mlonp,mlatp
  integer, pointer :: llon,llat

  ! named pointers for ease of use with external programs
  real(wp), pointer :: flagdirich      ! this will need to be converted to integer at some point because all datasets are floating point here
  real(wp), dimension(:,:), pointer :: E0xp,E0yp
  real(wp), dimension(:,:), pointer :: Vminx1p,Vmaxx1p
  real(wp), dimension(:), pointer :: Vminx2pslice,Vmaxx2pslice    !only slices because field lines (x1-dimension) are equipotentials
  real(wp), dimension(:), pointer :: Vminx3pslice,Vmaxx3pslice

  real(wp), dimension(:,:), pointer :: E0xiprev,E0yiprev
  real(wp), dimension(:,:), pointer :: Vminx1iprev,Vmaxx1iprev
  real(wp), dimension(:), pointer :: Vminx2isprev,Vmaxx2isprev    !only slices because field lines (x1-dimension) are equipotentials
  real(wp), dimension(:), pointer :: Vminx3isprev,Vmaxx3isprev

  real(wp), dimension(:,:), pointer :: E0xinext,E0yinext
  real(wp), dimension(:,:), pointer :: Vminx1inext,Vmaxx1inext
  real(wp), dimension(:), pointer :: Vminx2isnext,Vmaxx2isnext    !only slices because field lines (x1-dimension) are equipotentials
  real(wp), dimension(:), pointer :: Vminx3isnext,Vmaxx3isnext

  real(wp), dimension(:,:), pointer :: E0xinow,E0yinow
  real(wp), dimension(:,:), pointer :: Vminx1inow,Vmaxx1inow
  real(wp), dimension(:), pointer :: Vminx2isnow,Vmaxx2isnow    !only slices because field lines (x1-dimension) are equipotentials
  real(wp), dimension(:), pointer :: Vminx3isnow,Vmaxx3isnow
  contains
    ! overriding procedures
    procedure :: set_sizes=>set_sizes_efield

    ! deferred bindings
    procedure :: init=>init_efield
    procedure :: set_coordsi=>set_coordsi_efield
    procedure :: load_data=>load_data_efield
    procedure :: load_grid=>load_grid_efield
    procedure :: load_size=>load_size_efield
    final :: destructor
end type efielddata

contains
  !> need to override set_sizes so to account for fact that target interpolation is to a 2D grid, by default object will
  !    assume 3D and get the sizes from the simulation grid
  subroutine set_sizes_efield(self, &
                     l0D, &
                     l1Dax1,l1Dax2,l1Dax3, &
                     l2Dax23,l2Dax12,l2Dax13, &
                     l3D, &
                     x)
    class(efielddata), intent(inout) :: self
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

    ! coordinate axis sizes for interpolation states; noting that for a dipole we (internally) permute
    !   the interpolation sites in the same sense.
    !select type (x)
    !  class is (dipolemesh)
    !    print*, ' efielddata:  detected dipole mesh'
    !    self%lc1i=x%lx1;       ! note this dataset has 1D and 2D target interpolation grid
    !    self%lc2i=x%lx3all; self%lc3i=x%lx2all;
    !    ! dipolemesh mesh permuted ~alt,lat,lon, whereas inputdata organized lon,lat
    !  class default
        self%lc1i=x%lx1;       ! note this dataset has 1D and 2D target interpolation grid
        self%lc2i=x%lx2all; self%lc3i=x%lx3all;
    !end select

    ! check that the user is trying something sensible
    !if (self%lc1==1 .and. self%lc1i/=1 .or. self%lc2==1 .and. self%lc2i/=1 &
    !        .or. self%lc3==1 .and. self%lc3i/=1) then
    !  error stop 'inputdata:set_sizes() - singleton dimensions must be same for source and destination.'
    !end if

    ! flag sizes as assigned
    self%flagsizes=.true.
  end subroutine set_sizes_efield


  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_efield(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(efielddata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg                 ! gemini config type for optional params
    character(*), intent(in) :: sourcedir    ! directory for precipitation data input
    class(curvmesh), intent(in) :: x                    ! curvmesh object
    real(wp), intent(in) :: dtmodel,dtdata                      ! model time step and cadence for input data from config.nml
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='electric field boundary conditions'
    call self%set_name(strname)
    self%flagdoinput=cfg%flagE0file/=0

    ! read the simulation size from the source directory and allocate arrays
    allocate(self%lc1,self%lc2,self%lc3)      ! these are pointers
    self%llon=>self%lc2; self%llat=>self%lc3;
    call self%load_size()
    call self%set_sizes(1, &
                       0,2,2, &
                       4,0,0, &
                       0, &
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%mlonp=>self%coord2; self%mlatp=>self%coord3;
    call self%load_grid()

    ! set input data array pointers to faciliate easy to read input code
    self%flagdirich=>self%data0D(1)
    self%E0xp=>self%data2Dax23(:,:,1)
    self%E0yp=>self%data2Dax23(:,:,2)
    self%Vminx1p=>self%data2Dax23(:,:,3)
    self%Vmaxx1p=>self%data2Dax23(:,:,4)
    self%Vminx2pslice=>self%data1Dax3(:,1)
    self%Vmaxx2pslice=>self%data1Dax3(:,2)
    self%Vminx3pslice=>self%data1Dax2(:,1)
    self%Vmaxx3pslice=>self%data1Dax2(:,2)

    self%E0xiprev=>self%data2Dax23i(:,:,1,1)
    self%E0yiprev=>self%data2Dax23i(:,:,2,1)
    self%Vminx1iprev=>self%data2Dax23i(:,:,3,1)
    self%Vmaxx1iprev=>self%data2Dax23i(:,:,4,1)
    self%Vminx2isprev=>self%data1Dax3i(:,1,1)
    self%Vmaxx2isprev=>self%data1Dax3i(:,2,1)
    self%Vminx3isprev=>self%data1Dax2i(:,1,1)
    self%Vmaxx3isprev=>self%data1Dax2i(:,2,1)

    self%E0xinext=>self%data2Dax23i(:,:,1,2)
    self%E0yinext=>self%data2Dax23i(:,:,2,2)
    self%Vminx1inext=>self%data2Dax23i(:,:,3,2)
    self%Vmaxx1inext=>self%data2Dax23i(:,:,4,2)
    self%Vminx2isnext=>self%data1Dax3i(:,1,2)
    self%Vmaxx2isnext=>self%data1Dax3i(:,2,2)
    self%Vminx3isnext=>self%data1Dax2i(:,1,2)
    self%Vmaxx3isnext=>self%data1Dax2i(:,2,2)

    self%E0xinow=>self%data2Dax23inow(:,:,1)
    self%E0yinow=>self%data2Dax23inow(:,:,2)
    self%Vminx1inow=>self%data2Dax23inow(:,:,3)
    self%Vmaxx1inow=>self%data2Dax23inow(:,:,4)
    self%Vminx2isnow=>self%data1Dax3inow(:,1)
    self%Vmaxx2isnow=>self%data1Dax3inow(:,2)
    self%Vminx3isnow=>self%data1Dax2inow(:,1)
    self%Vmaxx3isnow=>self%data1Dax2inow(:,2)

    ! initialize the prev variables sensibly
    self%E0xiprev=0.0
    self%E0yiprev=0.0
    self%Vminx1iprev=0.0
    self%Vmaxx1iprev=0.0
    self%Vminx2isprev=0.0
    self%Vmaxx2isprev=0.0
    self%Vminx3isprev=0.0
    self%Vmaxx3isprev=0.0

    ! set to start time of simulation - set the first time update is called
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    print*, 'Preparing to prime input data arrays...'
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_efield


  !> get the input grid size from file, all workers will just call this sicne this is a one-time thing
  subroutine load_size_efield(self)
    class(efielddata), intent(inout) :: self

    ! basic error checking
    if (.not. self%flagsource) error stop 'efielddata:load_size_precip() - must define a source directory first'

    ! read sizes
    print '(/,A,/,A)', 'Electric field input:','--------------------'
    print '(A)', 'READ electric field size from: ' // self%sourcedir
    call get_simsize2(self%sourcedir // "/simsize.h5", llon=self%llon, llat=self%llat)

    print '(A,2I6)', 'Electric field size: llon,llat:  ',self%llon,self%llat
    if (self%llon < 1 .or. self%llat < 1) then
     print*, '  efielddata grid size must be strictly positive: ' //  self%sourcedir
     error stop
    end if

    ! set dim 1 size to null since inherently not used
    self%lc1=0

    ! flag to denote input data size is set
    self%flagdatasize=.true.
  end subroutine load_size_efield


  !> get the grid information from a file, all workers will just call this since one-time
  subroutine load_grid_efield(self)
    class(efielddata), intent(inout) :: self

    ! read grid data
    call get_grid2(self%sourcedir // "/simgrid.h5", self%mlonp, self%mlatp)

    print '(A,4F9.3)', 'Electric field mlon,mlat extent:  ',minval(self%mlonp(:)),maxval(self%mlonp(:)), &
                                                           minval(self%mlatp(:)),maxval(self%mlatp(:))
    if(.not. all(ieee_is_finite(self%mlonp))) error stop 'efielddata:loadgrid() - mlon must be finite'
    if(.not. all(ieee_is_finite(self%mlatp))) error stop 'efielddata:loadgrid() - mlat must be finite'
  end subroutine load_grid_efield


  subroutine set_coordsi_efield(self,cfg,x)
    class(efielddata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg     ! presently not used but possibly eventually?
    class(curvmesh), intent(in) :: x
    integer :: ix2,ix3,iflat,ix1ref,ix2ref,ix3ref,ix1offset,ix2offset,ix3offset

    ! source arrays in the grid object may have ghost cells; these are offsets for arrays that do not
    !   preserve lbound and ubound, e.g. array(:,1,1) and the like
    ix1offset=1-lbound(x%rall,1)
    ix2offset=1-lbound(x%rall,2)
    ix3offset=1-lbound(x%rall,3)

    !! reference locations for determining points onto which we are interpolating
    !! these are grid specific, not object specific...
    if (x%lx2all > 1 .and. x%lx3all>1) then ! 3D sim
      ix2ref = x%lx2all/2      !note integer division
      ix3ref = x%lx3all/2
    else if (x%lx2all==1 .and. x%lx3all>1) then
      ix2ref = 1
      ix3ref=x%lx3all/2
    else if (x%lx2all>1 .and. x%lx3all==1) then
      ix2ref=x%lx2all/2
      ix3ref=1
    else
      error stop 'Unable to orient boundary conditions for electric potential'
    endif

    !! by default the code uses 300km altitude as a reference location, using the center x2,x3 point
    !! These are the coordinates for inputs varying along axes 2,3
    ix1ref = minloc(abs(x%rall(:,ix2ref,ix3ref) - Re - 300e3_wp), dim=1)    ! includes ghost cells if x%rall has ghost cells
    ix1ref=ix1ref-ix1offset
    do ix3=1,lx3all
      do ix2=1,lx2all
        iflat=(ix3-1)*lx2all+ix2
        self%coord2iax23(iflat)=x%phiall(ix1ref,ix2,ix3)*180/pi
        self%coord3iax23(iflat)=90-x%thetaall(ix1ref,ix2,ix3)*180/pi
      end do
    end do
    if (debug) print '(A,4F7.2)', 'Grid has mlon,mlat range:  ',minval(self%coord2iax23),maxval(self%coord2iax23), &
                                     minval(self%coord3iax23),maxval(self%coord3iax23)
    if (debug) print *, 'Grid has size:  ',iflat

    if (self%flagdipmesh) then
      !! for electric field input data we also have some things that vary along axis 3 only
      do ix2=1,x%lx2all    ! note mangling ix2->ix3
        self%coord3iax3(ix2)=90-x%thetaall(ix1ref,ix2,1)*180/pi     ! default to ix2=1 side of the grid
      end do
      !! for BCs varing along axis 2 only
      do ix3=1,x%lx3all    ! note mangling ix3->ix2
        self%coord2iax2(ix3)=x%phiall(ix1ref,1,ix3)*180/pi          ! default to ix3=1 side of the grid
      end do
    else
      !! for electric field input data we also have some things that vary along axis 3 only
      do ix3=1,x%lx3all
        self%coord3iax3(ix3)=90-x%thetaall(ix1ref,1,ix3)*180/pi     ! default to ix2=1 side of the grid
      end do
      !! for BCs varing along axis 2 only
      do ix2=1,x%lx2all
        self%coord2iax2(ix2)=x%phiall(ix1ref,ix2,1)*180/pi          ! default to ix3=1 side of the grid
      end do
    end if

    !! mark coordinates as set
    self%flagcoordsi=.true.
  end subroutine set_coordsi_efield


  !> have root read in data from a file (only root manipulates this object)
  subroutine load_data_efield(self,t,dtmodel,ymdtmp,UTsectmp)
    class(efielddata), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
    integer :: flagdirich_int

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warning

    !! all workers should update the date
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt, ymdtmp, UTsectmp)

    !! all workers read data out of this file
    print*, '  date and time:  ',ymdtmp,UTsectmp
    print*, '  efield filename:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp)

    call get_Efield(date_filename(self%sourcedir, ymdtmp, UTsectmp) // ".h5", &
      flagdirich_int,self%E0xp,self%E0yp,self%Vminx1p,self%Vmaxx1p,&
      self%Vminx2pslice,self%Vmaxx2pslice,self%Vminx3pslice,self%Vmaxx3pslice)
    self%flagdirich=real(flagdirich_int,wp)

    if (debug) then
      print*, ' Solve type:  ', self%flagdirich
      print *, 'Min/max values for E0xp:  ',minval(self%E0xp),maxval(self%E0xp)
      print *, 'Min/max values for E0yp:  ',minval(self%E0yp),maxval(self%E0yp)
      print *, 'Min/max values for Vminx1p:  ',minval(self%Vminx1p),maxval(self%Vminx1p)
      print *, 'Min/max values for Vmaxx1p:  ',minval(self%Vmaxx1p),maxval(self%Vmaxx1p)
      print *, 'Min/max values for Vminx2pslice:  ',minval(self%Vminx2pslice),maxval(self%Vminx2pslice)
      print *, 'Min/max values for Vmaxx2pslice:  ',minval(self%Vmaxx2pslice),maxval(self%Vmaxx2pslice)
      print *, 'Min/max values for Vminx3pslice:  ',minval(self%Vminx3pslice),maxval(self%Vminx3pslice)
      print *, 'Min/max values for Vmaxx3pslice:  ',minval(self%Vmaxx3pslice),maxval(self%Vmaxx3pslice)
    endif

    if (.not. all(ieee_is_finite(self%E0xp))) error stop 'E0xp: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%E0yp))) error stop 'E0yp: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vminx1p))) error stop 'Vminx1p: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vmaxx1p))) error stop 'Vmaxx1p: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vminx2pslice))) error stop 'Vminx2pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vmaxx2pslice))) error stop 'Vmaxx2pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vminx3pslice))) error stop 'Vminx3pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(self%Vmaxx3pslice))) error stop 'Vmaxx3pslice: non-finite value(s)'
  end subroutine load_data_efield


  !> destructor needs to clear memory out
  subroutine destructor(self)
    type(efielddata), intent(inout) :: self

    call self%dissociate_pointers()
  end subroutine destructor
end module efielddataobj
