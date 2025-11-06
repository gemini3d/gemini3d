module neutraldataBGobj

! type extension for file-based neutral background atmospheric data from a profile

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use phys_consts, only: wp,debug,pi
use inputdataobj, only: inputdata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize3,get_grid3,get_neutralBG
use timeutils, only: dateinc,date_filename

implicit none (type, external)
private
public :: neutraldataBG

type, extends(inputdata) :: neutraldataBG
  ! coordinates for input precipitation data, and storage
  real(wp), dimension(:), pointer :: altp,glonp,glatp
  integer, pointer :: llon,llat,lalt
  real(wp), dimension(:,:,:,:), pointer :: natmp
  real(wp), dimension(:,:,:,:), pointer :: natmiprev,natminext,natminow

  ! work and target coordinates
  real(wp), dimension(:), pointer :: alti,gloni,glati

  ! projection factors needed to rotate input data onto grid
  real(wp), dimension(:,:,:), allocatable :: proj_ezp_e1,proj_ezp_e2,proj_ezp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_exp_e1,proj_exp_e2,proj_exp_e3

  contains
    ! deferred bindings
    procedure :: init=>init_neutralBG
    procedure :: set_coordsi=>set_coordsi_neutralBG
    procedure :: load_data=>load_data_neutralBG
    procedure :: load_grid=>load_grid_neutralBG
    procedure :: load_size=>load_size_neutralBG     ! load the size of the input data files

    ! overriding procedures
    !procedure :: update

    ! unique to this class
    !procedure :: rotate_winds

    ! final
    final :: destructor
end type neutraldataBG

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem) and prime everything
  !    so we are ready to call self%update()
  !  After this procedure is called all pointer aliases are set and can be used; internal to this procedure pay attention
  !  to ordering of when pointers are set with respect to when various type-bound procedures are called
  subroutine init_neutralBG(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutraldataBG), intent(inout) :: self
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
                       9, &     ! target data for neutralBG info is a 3D set of arrays
                       x )
    call self%init_storage()
    call self%set_cadence(dtdata)

    ! set local pointers grid pointers and assign input data grid
    self%altp=>self%coord1
    call self%load_grid()

    ! Set input data array pointers to faciliate easy to read input code; these may or may not be helpful to user
    self%natmp=>self%data3D(:,:,:,:)   !contiguous in memory since only one element in x2,3
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
    class(neutraldataBG), intent(inout) :: self
    integer :: ltmp    ! throwaway variable

    ! basic error checking
    if (.not. self%flagsource) error stop 'neutraldataBG:load_size_neutralBG() - must define a source directory first'

    ! read sizes
    call get_simsize3(self%sourcedir // "/simsize.h5", self%lalt, self%llon, self%llat)

    if (self%lalt < 1) then
     print*, '  neutralBG grid size must be strictly positive: ' //  self%sourcedir
     error stop
    end if

    ! flag to denote input data size is set
    self%flagdatasize=.true.
  end subroutine load_size_neutralBG


  !> get the grid information from a file, all workers will just call this since one-time
  subroutine load_grid_neutralBG(self)
    class(neutraldataBG), intent(inout) :: self

    ! read grid data
    call get_grid3(self%sourcedir // "/simgrid.h5", self%altp, self%glonp, self%glatp)

    if(.not. all(ieee_is_finite(self%altp))) error stop 'neutralBGBCs_fileinput: alt must be finite'
    if(.not. all(ieee_is_finite(self%glonp))) error stop 'neutralBGBCs_fileinput: glon must be finite'
    if(.not. all(ieee_is_finite(self%glatp))) error stop 'neutralBGBCs_fileinput: glat must be finite'
  end subroutine load_grid_neutralBG


  !> set target coordinates for interpolation sights
  subroutine set_coordsi_neutralBG(self,cfg,x)
    class(neutraldataBG), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg     ! presently not used but possibly eventually?
    class(curvmesh), intent(in) :: x
    integer :: ix1,ix2,ix3
    real(wp), dimension(:,:,:), allocatable :: altimat,glonimat,glatimat
    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3,3) :: ealt,eglat,eglon   
    real(wp) :: tmpsca
    real(wp), dimension(3) :: tmpvec, exprm, eyp, ezp

    ! aliases for target interpolation sites
    self%alti=>self%coord1i
    self%gloni=>self%coord2i
    self%glati=>self%coord3i

    allocate(altimat(1:x%lx1,1:x%lx2,1:x%lx3))
    allocate(glonimat,glatimat,mold=altimat)

    ! Target coordinates are 3D geographic
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          altimat(ix1,ix2,ix3)=x%alt(ix1,ix2,ix3)
          glonimat(ix1,ix2,ix3)=x%glon(ix1,ix2,ix3)
          glatimat(ix1,ix2,ix3)=x%glat(ix1,ix2,ix3)
        end do
      end do
    end do
    self%alti=pack(altimat,.true.)
    self%gloni=pack(glonimat,.true.)
    self%glati=pack(glatimat,.true.)
    deallocate(altimat,glonimat,glatimat)

    ! Storage of projections needed to rotate winds into geographic components (alt,glon,glat) -> (z,x,y)
    allocate(self%proj_ezp_e1(x%lx1,x%lx2,x%lx3),self%proj_ezp_e2(x%lx1,x%lx2,x%lx3),self%proj_ezp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_eyp_e1(x%lx1,x%lx2,x%lx3),self%proj_eyp_e2(x%lx1,x%lx2,x%lx3),self%proj_eyp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_exp_e1(x%lx1,x%lx2,x%lx3),self%proj_exp_e2(x%lx1,x%lx2,x%lx3),self%proj_exp_e3(x%lx1,x%lx2,x%lx3))

    call x%calc_unitvec_geo(ealt,eglon,eglat)   
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections as well)
          ezp=ealt(ix1,ix2,ix3,:)
          !ezp=x%er(ix1,ix2,ix3,:)

          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca

          ! we now need geographic unit vectors which we can get from our grid methods
          eyp=eglat(ix1,ix2,ix3,:)
          !eyp= -x%etheta(ix1,ix2,ix3,:)

          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca

          exprm=eglon(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)
          !exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)

          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e3(ix1,ix2,ix3)=tmpsca
        end do
      end do
    end do

    self%flagcoordsi=.true.
  end subroutine set_coordsi_neutralBG


  !> have all processes read in data from file to avoid any message passing
  subroutine load_data_neutralBG(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutraldataBG), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warnings

    !! all workers should update the date
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt, ymdtmp, UTsectmp)

    !!!!!!  read in solar neutral background data from file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! this read must be done repeatedly through simulation so have only root do file io
    print*, '  date and time:  ',ymdtmp,UTsectmp
    print*, '  neutralBG filename:  ',date_filename(self%sourcedir,ymdtmp,UTsectmp)

    ! read in the data for the "next" frame from file
!    call get_neutralBG(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", &
!            self%nOp,self%nN2p,self%nO2,self%nH,self%nN,self%vnx,self%vny,self%Tn)
    call get_neutralBG(date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5", &
            self%natmp(:,:,:,1),self%natmp(:,:,:,2),self%natmp(:,:,:,3),self%natmp(:,:,:,4), &
            self%natmp(:,:,:,5),self%natmp(:,:,:,7),self%natmp(:,:,:,8), &
            self%natmp(:,:,:,9) )
    self%natmp(:,:,:,6)=0._wp    ! vertical drift
    !print*, 'min/max data:  ',  minval(self%Iinfp),maxval(self%Iinfp)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine load_data_neutralBG


!    !> overriding procedure for updating neutral atmos (need additional rotation steps)
!  subroutine update(self,cfg,dtmodel,t,x,ymd,UTsec)
!    class(neutraldata3D_mpi), intent(inout) :: self
!    type(gemini_cfg), intent(in) :: cfg
!    real(wp), intent(in) :: dtmodel             ! need both model and input data time stepping
!    real(wp), intent(in) :: t                   ! simulation absoluate time for which perturabation is to be computed
!    class(curvmesh), intent(in) :: x            ! mesh object
!    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
!    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations
!
!    ! execute a basic update
!    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)
!
!    call self%rotate_winds()
!  end subroutine update
!
!
!  ! FIXME: the model already assumes the background data will be in geographic coordinates!!!
!  !> This subroutine takes winds stored in self%dvn?inow and applies a rotational transformation onto the
!  !      grid object for this simulation
!  subroutine rotate_winds(self)
!    class(neutraldata3D), intent(inout) :: self
!    integer :: ix1,ix2,ix3
!    real(wp) :: vnx,vny,vnz
!
!    ! do rotations one grid point at a time to cut down on temp storage needed
!    do ix3=1,self%lc3i
!      do ix2=1,self%lc2i
!        do ix1=1,self%lc1i
!          vnz=self%dvn1inow(ix1,ix2,ix3)    ! geographic altitude direction prior to rotation
!          vnx=self%dvn2inow(ix1,ix2,ix3)    ! geographic east (longitude) prior to rotation
!          vny=self%dvn3inow(ix1,ix2,ix3)    ! geograhpic north (latitude) prior to rotation
!          self%dvn1inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e1(ix1,ix2,ix3) + vnx*self%proj_exp_e1(ix1,ix2,ix3) + &
!                                        vny*self%proj_eyp_e1(ix1,ix2,ix3)
!          self%dvn2inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e2(ix1,ix2,ix3) + vnx*self%proj_exp_e2(ix1,ix2,ix3) + &
!                                        vny*self%proj_eyp_e2(ix1,ix2,ix3)
!          self%dvn3inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e3(ix1,ix2,ix3) + vnx*self%proj_exp_e3(ix1,ix2,ix3) + &
!                                        vny*self%proj_eyp_e3(ix1,ix2,ix3)
!        end do
!      end do
!    end do
!  end subroutine rotate_winds


  !> destructor needs to clear memory out
  subroutine destructor(self)
    type(neutraldataBG), intent(inout) :: self

    if (self%flagcoordsi) then
      ! in addition to the normal coordsi allocatables we also have projections for this extension
      deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
      deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
      deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
    end if

    call self%dissociate_pointers()
  end subroutine destructor
end module neutraldataBGobj
