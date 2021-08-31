module neutraldata3Dobj

use phys_consts, only: wp
use neutraldataobj, only: neutraldata
implicit none (type,external)

!> type definition for 3D neutral data
type, extends(neutraldata) :: neutraldata3D
  ! source data coordinate pointers
  real(wp), dimension(:), pointer :: xn,yn,zn
  integer, pointer :: lxn,lyn,lzn

  ! source data pointers
  real(wp), dimension(:,:,:), pointer ::

  ! projection factors needed to rotate input data onto grid
  real(wp), dimension(:,:,:), allocatable :: proj_ezp_e1,proj_ezp_e2,proj_ezp_e3    !these projections are used in the axisymmetric interpolation
  real(wp), dimension(:,:,:), allocatable :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3    !these are for Cartesian projections
  real(wp), dimension(:,:,:), allocatable :: proj_exp_e1,proj_exp_e2,proj_exp_e3 
  contains
    ! bindings for deferred procedures
    procedure :: init=>init_neu3D
    procedure :: load_data=>load_data_neu3D
    procedure :: load_grid=>load_grid_neu3D
    procedure :: load_size=>load_size_neu3D
    procedure :: set_coordsi=>set_coordsi_neu3D
end type neutraldata3D

contains
  !> create storage for arrays needed specifically for 3D neutral input calculations
  subroutine init_storage_neu3D(self)
    class(neutraldata3D), intent(inout) :: self

    ! geometric projections for vectors
    
  end subroutine init_storage_neu3D



  !> initialize storage for this type of neutral input data
  subroutine init_neu3D(self,sourcedir,dtneu)
    class(neutraldata3D), intent(inout) :: self
    character, dimension(:), intent(in) :: sourcedir
    integer :: lc1,lc2,lc3
    
    ! read in neutral grid size from file (io module)
    call getsimsize(sourcedir,lc1,lc2,lc3)

    ! set sizes, we have 7 arrays all 3D (irrespective of 2D vs. 3D neutral input)
    allocate(self%lc1,self%lc2,self%lc3)      ! these are pointers
    self%lzn=>lc1; self%lxn=>lc2; self%lyn=>lc3;
    call self%load_size()
    call self%set_sizes(lc1,lc2,lc3, &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data
             7)            ! number 3D data

    ! allocate space for arrays
    call self%init_storage()
    call self%init_stoarge_neu3D()     ! additional work array space for 3D neutral input
    call self%set_cadence(dtneu)

    ! alias coordinate arrays for internal calculations (as needed)
    self%zn=>self%coord1; self%xn=>self%coord2; self%yn=>self%coord3;
    call self%load_grid()

    ! set aliases to point to correct source data arrays
    self%dnO=>self%data3D(:,:,:,1)
    self%dnN2=>self%data3D(:,:,:,2)
    self%dnO2=>self%data3D(:,:,:,3)
    self%dvnz=>self%data3D(:,:,:,4)
    self%dvnx=>self%data3D(:,:,:,5)
    self%dvny=>self%data3D(:,:,:,6)
    self%dTn=>self%data3D(:,:,:,7)

    ! call to base class procedure to set pointers for prev,now,next
    call self%setptrs_grid()

    ! initialize previous data so we get a correct starting value
    self%dnOiprev=0
    self%dnN2iprev=0
    self%dnO2iprev=0
    self%dvn1iprev=0
    self%dvn2iprev=0
    self%dvn3iprev=0
    self%dTniprev=0

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_neu3D


  !> set coordinates for target interpolation points
  subroutine set_coordsi_neu3D(self,cfg,x)
    class(neutraldata3D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: meanyn
    real(wp) :: meanxn
    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
    real(wp) :: xp,yp
    real(wp), dimension(3) :: erhop,ezp,eyp,tmpvec,exprm
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid,ierr
    real(wp), dimension(x%lx1,x%lx2,x%lx3) :: zimat,rhoimat,yimat,ximat
    real(wp) :: maxzn
    real(wp), dimension(2) :: xnrange,ynrange
    integer, dimension(6) :: indices
    real(wp), dimension(:), pointer :: xi,yi,zi


    ! Alias target locations for clarity of code
    zi=>self%coord1i; xi=>self%coord2i; yi=>self%coord3i;


    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2 - cfg%sourcemlat*pi/180
    
    
    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
    if (mpi_cfg%myid==0) then
      print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
    end if
    zimat=x%alt     !vertical coordinate
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          !INTERPOLATION BASED ON GEOMAGNETIC COORDINATES
          theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle
          if (lx2/=1) then
            phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
          else
            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
          end if
    
    
          !!COMPUTE DISTANCES - ZZZ possibly superfluous for 3D case???
          !gammarads=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)     !this is actually cos(gamma)
          !if (gammarads>1) then     !handles weird precision issues in 2D
          !  gammarads=1
          !else if (gammarads<-1) then
          !  gammarads=-1
          !end if
          !gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
          !rhoimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''
          !! ZZZ end possibly superfluous block of code...
    
          !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
          theta3=theta2
          phi3=phi1
          gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
          if (gamma1 > 1) then     !handles weird precision issues in 2D
            gamma1 = 1
          else if (gamma1 < -1) then
            gamma1 = -1
          end if
          gamma1=acos(gamma1)
    
          gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
          if (gamma2 > 1) then     !handles weird precision issues in 2D
            gamma2= 1
          else if (gamma2 < -1) then
            gamma2= -1
          end if
          gamma2=acos(gamma2)
          xp=Re*gamma1
          yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement
    
    
          !COMPUTE COORDINATES FROM DISTANCES
          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
            yp= -yp            !do we want an abs here to be safe
          end if
          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
            xp= -xp
          end if
          !phip=atan2(yp,xp)
    
          ximat(ix1,ix2,ix3)=xp     !eastward distance
          yimat(ix1,ix2,ix3)=yp     !northward distance
    
    
          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
          ezp=x%er(ix1,ix2,ix3,:)
    
          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_ezp_e2(ix1,ix2,ix3)=tmpsca
    
          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_ezp_e1(ix1,ix2,ix3)=tmpsca
    
          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          proj_ezp_e3(ix1,ix2,ix3)=tmpsca
    
          eyp= -x%etheta(ix1,ix2,ix3,:)
    
          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_eyp_e1(ix1,ix2,ix3)=tmpsca
    
          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_eyp_e2(ix1,ix2,ix3)=tmpsca
    
          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_eyp_e3(ix1,ix2,ix3)=tmpsca
    
          exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)
    
          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_exp_e1(ix1,ix2,ix3)=tmpsca
    
          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_exp_e2(ix1,ix2,ix3)=tmpsca
    
          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          proj_exp_e3(ix1,ix2,ix3)=tmpsca
    
        end do
      end do
    end do
    
    
    !Assign values for flat lists of grid points
    zi=pack(zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
    yi=pack(yimat,.true.)
    xi=pack(ximat,.true.)

    ! FIXME: do we need to have the new grid code clear its unit vectors?
    if (mpi_cfg%myid==0) then
      print*, '...Clearing out unit vectors (after projections)...'
    end if
    !call clear_unitvecs(x)
    
    if(mpi_cfg%myid==0) then
      print*, 'Projection checking:  ',minval(proj_exp_e1),maxval(proj_exp_e1),minval(proj_exp_e2),maxval(proj_exp_e2), &
                                        minval(proj_exp_e3),maxval(proj_exp_e3)
    end if
  end subroutine set_coordsi_neu3D


  !> load the size information for the dataset we are dealing with
  subroutine load_size_neu3D(self)
    class(neutraldata3D), intent(inout) :: self

    ! basic error checking
    if (.not. self%flagsource) error stop 'neutraldata3D:load_size_neu3D() - must define a source directory first'

   !Establish the size of the grid based on input file and distribute to workers
if (mpi_cfg%myid==0) then    !root
  print '(A,/,A)', 'READ neutral size from:', cfg%sourcedir

  call get_simsize3(cfg%sourcedir, lx1=lxnall, lx2all=lynall, lx3all=lzn)

  print *, 'Neutral data has lx,ly,lz size:  ',lxnall,lynall,lzn,' with spacing dx,dy,dz',cfg%dxn,cfg%drhon,cfg%dzn
  if (lxnall < 1 .or. lynall < 1 .or. lzn < 1) then
    write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
    error stop 'neutral:gridproj_dneu3D: grid size must be strictly positive'
  endif

  !! 3D will not longer support storing fullgrid variables; wastes too much memory
  !allocate(dnOall(lzn,lxnall,lynall),dnN2all(lzn,lxnall,lynall),dnO2all(lzn,lxnall,lynall),dvnrhoall(lzn,lxnall,lynall), &
  !            dvnzall(lzn,lxnall,lynall),dvnxall(lzn,lxnall,lynall),dTnall(lzn,lxnall,lynall))    !ZZZ - note that these might be deallocated after each read to clean up memory management a bit...

  !> alert workers of the grid size in zn
  do iid=1,mpi_cfg%lid-1
    call mpi_send(lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
  end do
else
  call mpi_recv(lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if 

  end subroutine load_size_neu3D


  !> unlike other inputs the neutral input grid is computed from numnber of points and a step size
  subroutine load_grid_neu3D



!! ROOT
  !calculate the z grid (same for all) and distribute to workers so we can figure out their x-y slabs
  print*, '...creating vertical grid and sending to workers...'
  zn=[ ((real(izn, wp)-1)*cfg%dzn, izn=1,lzn) ]    !root calculates and distributes but this is the same for all workers - assmes that the max neutral grid extent in altitude is always less than the plasma grid (should almost always be true)
  maxzn=maxval(zn)

    call mpi_send(zn,lzn,mpi_realprec,iid,tag%zn,MPI_COMM_WORLD,ierr)





  end subroutine
end module neutraldata3Dobj
