module inputdataobj

use phys_consts, only : wp
use meshobj, only : curvmesh
use interpolation, only : interp1,interp2,interp3
use timeutils, only : dateinc, date_filename, find_lastdate

implicit none (type, external)
public


!> this is a generic class for an data object being input into the model and interpolated in space and time
type, abstract :: inputdata
  !! a name for our dataset
  character, dimension(:), allocatable :: dataname='DEFAULT'

  !! here we store data that have already been received but not yet interpolated
  real(wp), dimension(:), pointer :: coord1,coord2,coord3     ! coordinates for the source data (interpolant coords)
  integer :: lc1,lc2,lc3                                      ! dataset length along the 3 coordinate axes
  real(wp), dimension(:), pointer :: data0D
  real(wp), dimension(:,:), pointer :: data1Dax1,data1Dax2,data1Dax3
  real(wp), dimension(:,:), pointer :: data2Dax23,data2Dax12,data2Dax13
  real(wp), dimension(:,:,:), pointer :: data3D

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
  integer :: lparm3D
  real(wp), dimension(:,:,:), pointer :: coord1i,coord2i,coord3i     ! coordinates of the interpolation sites
  integer :: lc1i,lc2i,lc3i                                          ! dataset length along the 3 coordinate axes

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
    !! top-level user-intended
    procedure :: update                          ! check to see if new file needs to be read and read accordingly (will need to call deferred loaddata)
    procedure(makeproc), deferred :: make        ! get up object for first time step:  call init_storage, call prime_data, set data cadence

    !! internal/fine-grained control
    procedure :: set_sizes             ! initiate sizes for coordinate axes and number of datasets of different dimensionality
    procedure :: set_coords            ! fill interpolant coordinate arrays
    procedure :: set_cadence           ! fill dt variable for this instance of input
    procedure :: set_name              ! assign a character string name to our dataset
    procedure :: init_storage          ! wrapper routine to set up arrays once sizes are known/set
    procedure :: spaceinterp           ! interpolate spatially
    procedure :: timeinterp            ! interpolate in time based on data presently loaded into spatial arrays
    procedure :: dissociate_pointers   ! clear out memory and reset and allocation status flags
    procedure :: prime_data            ! load data buffers so that the object is ready for the first time step

    !! internal, data kind specific
    procedure(coordisetproc), deferred :: set_coordsi       ! use grid data to compute coordinates of the interpolation sites
    procedure(loadproc), deferred :: load_data             ! read data from file (possibly one array at a time) and spatially interpolate and store it in the appropriate arrays
end type inputdata


!> interfaces for deferred procedures
abstract interface
  subroutine makeproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine makeproc
  subroutine coordsetproc(self,cfg,x)
    import inputdata
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in)  :: x
  end subroutine coordisetproc
  subroutine loadproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine loadproc
end interface

contains
  !> "prime" data at the beginning of the simulation so that proper inputs can be derived/interpolated for the first time step
  !     Note that we need to separate any activity that isn't directly related to input data (e.g. background states, etc.) from 
  !     this routine so that it purely acts on properties of the inputdata class/type
  subroutine prime_data(self,cfg,x,dtdata)
    class(inputdata), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg      ! need start date, etc. of the simulation to prime
    class(curvmesh), intent(in) :: x         ! must have the grid to call update
    real(wp), intent(in) :: dtdata           ! need time step for *this* general data input
    integer, dimension(3) :: ymdprev,ymdnext   !dates for interpolated data
    real(wp) :: UTsecprev,UTsecnext

    !! find the last input data preceding the milestone/initial condition that we start with
    call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,dtdata,ymdtmp,UTsectmp)
  
    !! Loads the neutral input file corresponding to the "first" time step of the simulation to prevent the first interpolant
    !  from being zero and causing issues with restart simulations.  I.e. make sure the neutral buffers are primed for restart
    !  This requires us to load file input twice, once corresponding to the initial frame and once for the "first, next" frame.
    self%t(1)=UTsectmp-UTsec-2*dtdata
    self%t(2)=tprev+dtdata
    !if (mpi_cfg%myid==0) print*, '!!!Attempting initial load of neutral dynamics files!!!' // &
    !                         ' This is a workaround to insure compatibility with restarts...',ymdtmp,UTsectmp
    !! We essentially are loading up the data corresponding to halfway betwween -dtneu and t0 (zero).  This will load
    !   two time levels back so when tprev is incremented twice it will be the true tprev corresponding to first time step
    call self%update(cfg,dt,dtdata,self%t(2)+dtdata/2,ymdtmp,UTsectmp-dtdata, &
                          x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)  !abs time arg to be < 0
    ! FIXME: need to get rid of background arrays, etc. and have the neutral_update part done outside of this type/class
  
    !if (mpi_cfg%myid==0) print*, 'Now loading initial next file for neutral perturbations...'
    !! Now compute perturbations for the present time (zero), this moves the primed variables in next into prev and then
    !  loads up a current state so that we get a proper interpolation for the first time step.
    call self%update(cfg,dt,dtdata,0._wp,ymdtmp,UTsectmp,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)    !t-dt so we land exactly on start time
    ! FIXME: need to get rid of background arrays, etc. and have the neutral_update part done outside of this type/class
  end subroutine prime_data


  !> wrapper routine to execute various steps needed to update input data to present time
  subroutine update(self,t,dtmodel,sourcedir,x,ymd,UTsec)
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    character, dimension(:), intent(in) :: sourcedir
    class(curvmesh), intent(inout) :: x
    
    integer :: ix1,ix2,ix3,iid!,irhon,izn
    integer, dimension(3) :: ymdtmp
    real(wp) :: UTsectmp
    
    
    !! see if we need to load new data into the buffer; negative time means that we need to load the first frame
    if (t+dtmodel/2 >= t(2) .or. t < 0) then       
      !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
      if (.not. allocated(zn)) then     !means this is the first time we've tried to load neutral simulation data, should we check for a previous neutral file to load??? or just assume everything starts at zero?  This needs to somehow check for an existing file under certain conditiosn, maybe if it==1???  Actually we don't even need that we can just check that the neutral grid is allocated (or not)
        !initialize dates
        self%ymdref(:,1)=ymd
        self%UTsecref(1)=UTsec
        self%ymdref(:,2)=ymdref(:,1)
        self%UTsec(2)=UTsec(1)
    
        !Create a neutral grid, do some allocations and projections
        call self%set_coordi(cfg,x)    ! FIXME: how to differentiate between 2D cart vs. axisymmetric?  Need separate objects, extensions of an extension probably...
      end if
    
      !Read in neutral data from a file
      call self%load_data(t,dt,sourcedir,ymdtmp,UTsectmp)   !FIXME: need to differentiate between different types of 2D grids e.g. via subtype extensions perhaps?
    
      !Spatial interpolation for the frame we just read in
      if (mpi_cfg%myid==0 .and. debug) then
        print *, 'Spatial interpolation and rotations (if necessary) for dataset:  ',self%dataname,' for date:  ',ymdtmp,' ',UTsectmp
      end if
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
  end subroutine update


  !> use data stored in input arrays to interpolate onto grid sites for "next" dataset
  subroutine spaceinterp(self)
    class(inputdata),intent(inout) :: self
    integer :: iparm
    real(wp), dimension(:), allocatable :: tempdata

    !> 1D arrays varying along the 1-axis
    if (associated(data1Dax1i)) then
      allocate(tempdata(self%lc1i))
      do iparm=1,l1Dax1
        tempdata(:)=interp1(coord1,data1Dax1(:,iparm),coord1i)
        self%data1Dax1i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 2-axis
    if (associated(data1Dax2i)) then
      allocate(tempdata(self%lc2i))
      do iparm=1,l1Dax2
        tempdata(:)=interp1(coord2,data1Dax2(:,iparm),coord2i)
        self%data1Dax2i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 1D arrays varying along the 3-axis
    if (associated(data1Dax3i)) then
      allocate(tempdata(self%lc3i))
      do iparm=1,l1Dax3
        tempdata(:)=interp1(coord3,data1Dax3(:,iparm),coord3i)
        self%data1Dax3i(:,iparm,2)=tempdata(:)
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 2,3 axes
    if (associated(data2Dax23i)) then
      allocate(tempdata(self%lc2i*self%lc3i))
      do iparm=1,l2D23
        tempdata(:)=interp2(coord2,coord3,data2Dax23(:,:,iparm),coord2i,coord3i)
        self%data2Dax23i(:,:,iparm,2)=reshape(tempdata,[lc2i,lc3i])
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 1,2 axes
    if (associated(data2Dax12i)) then
      allocate(tempdata(self%lc1i*self%lc2i))
      do iparm=1,l2D12
        tempdata(:)=interp2(coord1,coord2,data2Dax12(:,:,iparm),coord1i,coord2i)
        self%data2Dax12i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i])
      end do
      deallocate(tempdata)
    end if

    !> 2D arrays varying along the 1,3 axes
    if (associated(data2Dax13i)) then
      allocate(tempdata(self%lc1i*self%lc3i))
      do iparm=1,l2D13
        tempdata(:)=interp2(coord1,coord3,data2Dax13(:,:,iparm),coord1i,coord3i)
        self%data2Dax13i(:,:,iparm,2)=reshape(tempdata,[lc1i,lc3i])
      end do
      deallocate(tempdata)
    end if

    !> 3D arrays varying along all axes (obv.)
    if (associated(data3Di)) then
      allocate(tempdata(self%lc1i*self%lc2i*self%lc3i))
      do iparm=1,l3D
        tempdata(:)=interp3(coord1,coord2,coord3,data3D(:,:,:,iparm),coord1i,coord2i,coord3i)
        self%data3Di(:,:,:,iparm,2)=reshape(tempdata,[lc1i,lc2i,lc3i])
      end do
      deallocate(tempdata)
    end if 
  end subroutine spaceinterp


  !> interpolate data in time (requires first loading and interpolating in space)
  subroutine timeinterp(self,t,dt)
    class(inputdata), intent(inout) :: self
    real(wp), intent(in) :: t
    real(wp) :: slope
    integer :: ic1,ic2,ic3

    !> interpolate scalars in time
    if (associated(data0Di)) then
      slope=(self%data0Di(iparm,2)-self%data0Di(iparm,1))/(self%tref(2)-self%tref(1))
      self%data0Dinow(ic1,iparm)=self%data0Di(iparm,1)+slope*(t+dt/2 -self%tref(1))
    end if

    !> 1D arrays varying along the 1-axis
    if (associated(data1Dax1i)) then
      do iparm=1,l1Dax1
        do ic1=1,lc1i
          slope=(self%data1Dax1i(ic1,iparm,2)-self%data1Dax1i(ic1,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax1inow(ic1,iparm)=self%data1Dax1i(ic1,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 1D arrays varying along the 2-axis
    if (associated(data1Dax2i)) then
      do iparm=1,l1Dax2
        do ic2=1,lc2i
          slope=(self%data1Dax2i(ic2,iparm,2)-self%data1Dax2i(ic2,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax2inow(ic2,iparm)=self%data1Dax2i(ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 1D arrays varying along the 3-axis
    if (associated(data1Dax3i)) then
      do iparm=1,l1Dax3
        do ic3=1,lc3i
          slope=(self%data1Dax3i(ic3,iparm,2)-self%data1Dax3i(ic3,iparm,1))/(self%tref(2)-self%tref(1))
          self%data1Dax3inow(ic3,iparm)=self%data1Dax3i(ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
        end do
      end do
    end if

    !> 2D arrays varying along the 2,3 axes
    if (associated(data2Dax23i)) then
      do iparm=1,l2D23
        do ic3=1,lc3i
          do ic2=1,lc2i
            slope=(self%data2Dax23i(ic2,ic3,iparm,2)-self%data2Dax23i(ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax23inow(ic2,ic3,iparm)=self%data2Dax23i(ic2,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 2D arrays varying along the 1,2 axes
    if (associated(data2Dax12i)) then
      do iparm=1,l2D12
        do ic2=1,lc2i
          do ic1=1,lc1i
            slope=(self%data2Dax12i(ic1,ic2,iparm,2)-self%data2Dax12i(ic1,ic2,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax12inow(ic1,ic2,iparm)=self%data2Dax12i(ic1,ic2,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 2D arrays varying along the 1,3 axes
    if (associated(data2Dax13i)) then
      do iparm=1,l2D13
        do ic3=1,lc3i
          do ic1=1,lc1i
            slope=(self%data2Dax13i(ic1,ic3,iparm,2)-self%data2Dax13i(ic1,ic3,iparm,1))/(self%tref(2)-self%tref(1))
            self%data2Dax13inow(ic1,ic3,iparm)=self%data2Dax13i(ic1,ic3,iparm,1)+slope*(t+dt/2 -self%tref(1))
          end do
        end do
      end do
    end if

    !> 3D arrays varying along all axes (obv.)
    if (associated(data3Di)) then
      do iparm=1,l3D
        do ic3=1,lc3i
          do ic2=1,lc2i
            do ic1=1,lc1i
              slope=(self%data3Di(ic1,ic2,ic3,iparm,2)-self%data3Di(ic1,ic2,ic3,iparm,1))/(self%tref(2)-self%tref(1))
              self%data3Dinow(ic1,ic2,ic3,iparm)=self%data3Di(ic1,ic2,ic3,iparm,1)+slope*(t+dt/2-self%tref(1))
            end do
          end do
        end do
      end do
    end if 

    !> update current time in object
    self%tnow=t+dt/2
  end subroutine timeinterp
end module inputdataobj
