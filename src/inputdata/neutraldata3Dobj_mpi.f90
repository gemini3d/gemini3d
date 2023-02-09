module neutraldata3Dobj_mpi

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize3,get_simsize2,get_grid2,get_precip
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi
use timeutils, only: dateinc,date_filename
use h5fortran, only: hdf5_file
use grid, only: gridflag

use mpi_f08, only : mpi_send,mpi_recv,mpi_integer,mpi_comm_world,mpi_status_ignore

implicit none (type, external)
private
public :: neutraldata3D_mpi


!> type definition for 3D neutral data that require mpi communication for loading data
type, abstract, extends(neutraldata3D) :: neutraldata3D_mpi
  ! mpi-related information on subgrid extents and indices, only used on the root process; otherwise ignored
  real(wp), dimension(:,:), allocatable :: extents    ! min/max x,y,z of each worker
  integer, dimension(:,:), allocatable :: indx        ! indices for each workers' pieces of the neutral data
  integer, dimension(:,:), allocatable :: slabsizes
  contains
    ! replacement for gridsize and gridload
    procedure :: load_sizeandgrid_neu3D

    ! overriding procedures
    procedure :: update

    ! bindings for deferred procedures
    procedure :: init=>init_neu3D
    procedure :: load_data=>load_data_neu3D
    procedure :: load_grid=>load_grid_neu3D    ! stub, does nothing since the grid information is handled differently for different
                                               ! extensions, viz. some will set this data/information in an different manner than the 
                                               ! usual inputdata objects
    procedure :: load_size=>load_size_neu3D    ! stub, child must override if used or keep no-op
end type neutraldata3D_mpi


!> interfaces for submodule "utility" procedures
interface ! neuslab.f90
  module subroutine slabrange(maxzn,ximat,yimat,zimat,sourcemlat,xnrange,ynrange,gridflag)
    real(wp), intent(in) :: maxzn
    real(wp), dimension(:,:,:), intent(in) :: ximat,yimat,zimat
    real(wp), intent(in) :: sourcemlat
    real(wp), dimension(2), intent(out) :: xnrange,ynrange     !for min and max
    integer, intent(in) :: gridflag
  end subroutine slabrange
  module subroutine  range2inds(ranges,zn,xnall,ynall,indices)
    real(wp), dimension(6), intent(in) :: ranges
    real(wp), dimension(:), intent(in) :: zn,xnall,ynall
    integer, dimension(6), intent(out) :: indices
  end subroutine range2inds
  module subroutine dneu_root2workers(paramall,tag,slabsizes,indx,param)
    real(wp), dimension(:,:,:), intent(in) :: paramall
    integer, intent(in) :: tag
    integer, dimension(0:,:), intent(in) :: slabsizes
    integer, dimension(0:,:), intent(in) :: indx
    real(wp), dimension(:,:,:), intent(inout) :: param
  end subroutine dneu_root2workers
  module subroutine dneu_workers_from_root(tag,param)
    integer, intent(in) :: tag
    real(wp), dimension(:,:,:), intent(inout) :: param
  end subroutine dneu_workers_from_root
end interface

contains
  !> initialize storage for this type of neutral input data
  subroutine init_neu3D(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutraldata3D_mpi), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel,dtdata
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation
    character(:), allocatable :: strname    ! allow auto-allocate for strings

    ! force 3D interpolation regardless of working subarray size
    self%flagforcenative=.true.

    ! tell our object where its data are and give the dataset a name
    call self%set_source(sourcedir)
    strname='neutral perturbations (3D)'
    call self%set_name(strname)
    call self%set_cadence(dtdata)
    self%flagdoinput=cfg%flagdneu/=0

    ! set sizes, we have 7 arrays all 3D (irrespective of 2D vs. 3D neutral input).  for 3D neutral input
    !    the situation is more complicated that for other datasets because you cannot compute the number of
    !    source grid points for each worker until you have root compute the entire grid and slice everything up
    allocate(self%lc1,self%lc2,self%lc3)                                     ! these are pointers, even though scalar
    self%lzn=>self%lc1; self%lxn=>self%lc2; self%lyn=>self%lc3;              ! these referenced while reading size and grid data
    call self%set_coordsi(cfg,x)                   ! since this preceeds init_storage it must do the work of allocating some spaces
    call self%load_sizeandgrid_neu3D(cfg)          ! cfg needed to form source neutral grid
    call self%set_sizes( &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data
             7, &          ! number 3D datasets
             x)          ! The main purpose of this is to set the number of 3D datasets (other params already set)

    ! allocate space for arrays, note for neutrals some of this has already happened so there is an overloaded procedure
    call self%init_storage()

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

    ! set to start time of simulation - not needed since assigned by update on first call.  FIXME: a bit messy
    !self%ymdref(:,1)=cfg%ymd0; self%ymdref(:,2)=cfg%ymd0;
    !self%UTsecref(1)=cfg%UTsec0; self%UTsecref(2)=cfg%UTsec0;

    ! prime input data
    call self%prime_data(cfg,x,dtmodel,ymd,UTsec)
  end subroutine init_neu3D



  !> load source data size and grid information and communicate to worker processes.
  !    Note that this routine will allocate sizes for source coordinates grids in constrast
  !    with other inputdata type extensions which have separate load_size, allocate, and
  !    load_grid procedures.
  subroutine load_sizeandgrid_neu3D(self,cfg)
    class(neutraldata3D_mpi), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:), allocatable :: xn,yn             ! for root to break off pieces of the entire grid array
    integer :: izn,iid                                  ! local copies for root, eventually these need to be stored in object
    real(wp) :: maxzn
    real(wp), dimension(2) :: xnrange,ynrange                ! these eventually get stored in extents
    integer, dimension(6) :: indices                         ! these eventually get stored in indx
    integer :: ixn,iyn
    integer :: lxn,lyn
    real(wp) :: meanxn,meanyn

    if (mpi_cfg%myid==0) then    !root must establish the size of the grid based on input file and distribute to workers
      print '(A,/,A)', 'READ neutral size from:', self%sourcedir
      call get_simsize3(self%sourcedir // "/simsize.h5", lx1=self%lxnall, lx2all=self%lynall, lx3all=self%lzn)
      print *, 'Neutral data has lx,ly,lz size:  ',self%lxnall,self%lynall,self%lzn, &
                   ' with spacing dx,dy,dz',cfg%dxn,cfg%drhon,cfg%dzn
      if (self%lxnall < 1 .or. self%lynall < 1 .or. self%lzn < 1) then
        write(stderr,*) 'ERROR: reading ' // self%sourcedir
        error stop 'neutral:gridproj_dneu3D: grid size must be strictly positive'
      endif

      ! allocate space for target coordinate and bind alias
      allocate(self%coord1(self%lzn))
      self%zn=>self%coord1
      allocate(self%xnall(self%lxnall))
      allocate(self%ynall(self%lynall))

      !calculate the z grid (same for all) and distribute to workers so we can figure out their x-y slabs
      print*, '...creating vertical grid and sending to workers...'
      self%zn=[ ((real(izn, wp)-1)*cfg%dzn, izn=1,self%lzn) ]    !root calculates and distributes but this is the same for all workers - assmes that the max neutral grid extent in altitude is always less than the plasma grid (should almost always be true)
      maxzn=maxval(self%zn)
      do iid=1,mpi_cfg%lid-1
        call mpi_send(self%lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD)
        call mpi_send(self%zn,self%lzn,mpi_realprec,iid,tag%zn,MPI_COMM_WORLD)
      end do

      !Define a neutral grid (input data) x,y extent by assuming that the spacing is constant
      self%ynall=[ ((real(iyn, wp)-1)*cfg%drhon, iyn=1,self%lynall) ]
      meanyn=sum(self%ynall,1)/size(self%ynall,1)
      self%ynall=self%ynall-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
      self%xnall=[ ((real(ixn, wp)-1)*cfg%dxn, ixn=1,self%lxnall) ]
      meanxn=sum(self%xnall,1)/size(self%xnall,1)
      self%xnall=self%xnall-meanxn     !the neutral grid should be centered on zero for a cartesian interpolation
      print *, 'Created full neutral grid with y,z extent:',minval(self%xnall),maxval(self%xnall),minval(self%ynall), &
                    maxval(self%ynall),minval(self%zn),maxval(self%zn)

      ! calculate the extents of root grid using max altitude specified for the neutral grid
      call slabrange(maxzn,self%ximat,self%yimat,self%zimat,cfg%sourcemlat,xnrange,ynrange,gridflag)
      allocate(self%extents(0:mpi_cfg%lid-1,6),self%indx(0:mpi_cfg%lid-1,6),self%slabsizes(0:mpi_cfg%lid-1,2))
      self%extents(0,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]

      !receive extents of each of the other workers: extents(mpi_cfg%lid,6)
      print*, 'Receiving xn and yn ranges from workers...'
      do iid=1,mpi_cfg%lid-1
        call mpi_recv(xnrange,2,mpi_realprec,iid,tag%xnrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
        call mpi_recv(ynrange,2,mpi_realprec,iid,tag%ynrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
        self%extents(iid,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]     !need to store values as xnrange overwritten for each worker
        print*, 'Subgrid extents:  ',iid,self%extents(iid,:)
      end do

      !find index into into neutral arrays for each worker:  indx(mpi_cfg%lid,6)
      print*, 'Root grid check:  ',self%ynall(1),self%ynall(self%lynall)
      print*, 'Converting ranges to indices...'
      do iid=0,mpi_cfg%lid-1
        call range2inds(self%extents(iid,1:6),self%zn,self%xnall,self%ynall,indices)
        self%indx(iid,1:6)=indices
        print*, 'Subgrid indices',iid,self%indx(iid,:)
      end do

      !send each worker the sizes for their particular chunk (all different) and send worker that grid chunk
      print*,'Sending sizes and xn,yn subgrids to workers...'
      do iid=1,mpi_cfg%lid-1
        lxn=self%indx(iid,4)-self%indx(iid,3)+1
        lyn=self%indx(iid,6)-self%indx(iid,5)+1
        self%slabsizes(iid,1:2)=[lxn,lyn]
        call mpi_send(lyn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD)
        call mpi_send(lxn,1,MPI_INTEGER,iid,tag%lx,MPI_COMM_WORLD)
        allocate(xn(lxn),yn(lyn))
        xn=self%xnall(self%indx(iid,3):self%indx(iid,4))
        yn=self%ynall(self%indx(iid,5):self%indx(iid,6))
        call mpi_send(xn,lxn,mpi_realprec,iid,tag%xn,MPI_COMM_WORLD)
        call mpi_send(yn,lyn,mpi_realprec,iid,tag%yn,MPI_COMM_WORLD)
        deallocate(xn,yn)
      end do

      !have root store its part to the full neutral grid
      print*, 'Root is picking out its own subgrid...'
      self%lxn=self%indx(0,4)-self%indx(0,3)+1
      self%lyn=self%indx(0,6)-self%indx(0,5)+1
      self%slabsizes(0,1:2)=[self%lxn,self%lyn]

      ! allocate space and bind alias
      allocate(self%coord2(self%lxn),self%coord3(self%lyn))
      self%xn=>self%coord2; self%yn=>self%coord3;        ! input data coordinates

      ! store source coordinates
      self%xn=self%xnall(self%indx(0,3):self%indx(0,4))
      self%yn=self%ynall(self%indx(0,5):self%indx(0,6))
    else                 !workers
      !get the z-grid from root so we know what the max altitude we have to deal with will be
      call mpi_recv(self%lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE)

      ! allocate space for target coordinate and bind alias
      allocate(self%coord1(self%lzn))
      self%zn=>self%coord1

      ! receive data from root
      call mpi_recv(self%zn,self%lzn,mpi_realprec,0,tag%zn,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      maxzn=maxval(self%zn)

      !calculate the extent of my grid
      call slabrange(maxzn,self%ximat,self%yimat,self%zimat,cfg%sourcemlat,xnrange,ynrange,gridflag)

      !send ranges to root
      call mpi_send(xnrange,2,mpi_realprec,0,tag%xnrange,MPI_COMM_WORLD)
      call mpi_send(ynrange,2,mpi_realprec,0,tag%ynrange,MPI_COMM_WORLD)

      !receive my sizes from root, allocate then receive my pieces of the grid
      call mpi_recv(self%lxn,1,MPI_INTEGER,0,tag%lx,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call mpi_recv(self%lyn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE)

      ! at this point we can allocate space for the source coordinates and bind aliases as needed
      allocate(self%coord2(self%lxn),self%coord3(self%lyn))
      self%xn=>self%coord2; self%yn=>self%coord3;        ! input data coordinates

      ! recieve data from root
      call mpi_recv(self%xn,self%lxn,mpi_realprec,0,tag%xn,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call mpi_recv(self%yn,self%lyn,mpi_realprec,0,tag%yn,MPI_COMM_WORLD,MPI_STATUS_IGNORE)
    end if

    self%flagdatasize=.true.
  end subroutine load_sizeandgrid_neu3D


!  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
!  subroutine set_coordsi_neu3D(self,cfg,x)
!    class(neutraldata3D), intent(inout) :: self
!    type(gemini_cfg), intent(in) :: cfg
!    class(curvmesh), intent(in) :: x
!    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
!    real(wp) :: xp,yp
!    real(wp), dimension(3) :: ezp,eyp,tmpvec,exprm
!    real(wp) :: tmpsca
!    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid
!
!
!    ! Space for coordinate sites and projections in neutraldata3D object
!    allocate(self%coord1i(x%lx1*x%lx2*x%lx3),self%coord2i(x%lx1*x%lx2*x%lx3),self%coord3i(x%lx1*x%lx2*x%lx3))
!    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i;     ! coordinates of interpolation sites
!    allocate(self%ximat(x%lx1,x%lx2,x%lx3),self%yimat(x%lx1,x%lx2,x%lx3),self%zimat(x%lx1,x%lx2,x%lx3))
!    allocate(self%proj_ezp_e1(x%lx1,x%lx2,x%lx3),self%proj_ezp_e2(x%lx1,x%lx2,x%lx3),self%proj_ezp_e3(x%lx1,x%lx2,x%lx3))
!    allocate(self%proj_eyp_e1(x%lx1,x%lx2,x%lx3),self%proj_eyp_e2(x%lx1,x%lx2,x%lx3),self%proj_eyp_e3(x%lx1,x%lx2,x%lx3))
!    allocate(self%proj_exp_e1(x%lx1,x%lx2,x%lx3),self%proj_exp_e2(x%lx1,x%lx2,x%lx3),self%proj_exp_e3(x%lx1,x%lx2,x%lx3))
!
!    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
!    phi1=cfg%sourcemlon*pi/180
!    theta1=pi/2 - cfg%sourcemlat*pi/180
!
!    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
!    if (mpi_cfg%myid==0) then
!      print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
!    end if
!
!    self%zimat=x%alt     !vertical coordinate is just altitude array already stored in grid object
!    do ix3=1,x%lx3
!      do ix2=1,x%lx2
!        do ix1=1,x%lx1
!          ! interpolation based on geomag
!          theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle
!
!          !print*, ' center NS set',shape(self%zi),shape(self%zimat),theta2,x%theta(ix1,ix2,ix3)
!
!          if (x%lx2/=1) then
!            phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
!          else
!            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
!          end if
!
!          !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
!          theta3=theta2
!          phi3=phi1
!          gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
!          if (gamma1 > 1) then     !handles weird precision issues in 2D
!            gamma1 = 1
!          else if (gamma1 < -1) then
!            gamma1 = -1
!          end if
!          gamma1=acos(gamma1)
!
!          gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
!          if (gamma2 > 1) then     !handles weird precision issues in 2D
!            gamma2= 1
!          else if (gamma2 < -1) then
!            gamma2= -1
!          end if
!          gamma2=acos(gamma2)
!          xp=Re*gamma1
!          yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement
!
!          ! coordinates from distances
!          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
!            yp= -yp            !do we want an abs here to be safe
!          end if
!          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
!            xp= -xp
!          end if
!
!          self%ximat(ix1,ix2,ix3)=xp     !eastward distance
!          self%yimat(ix1,ix2,ix3)=yp     !northward distance
!
!          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
!          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
!          ezp=x%er(ix1,ix2,ix3,:)
!
!          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
!          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca
!
!          eyp= -x%etheta(ix1,ix2,ix3,:)
!
!          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca
!
!          exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)
!
!          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_exp_e1(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_exp_e2(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_exp_e3(ix1,ix2,ix3)=tmpsca
!        end do
!      end do
!    end do
!
!    !Assign values for flat lists of grid points
!    if (mpi_cfg%myid==0) then
!      print*, '...Packing interpolation target points...'
!    end if
!    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
!    self%yi=pack(self%yimat,.true.)
!    self%xi=pack(self%ximat,.true.)
!
!    ! FIXME: do we need to have the new grid code clear its unit vectors?  Or maybe this isn't a huge waste of memory???
!    if (mpi_cfg%myid==0) then
!      print*, '...Clearing out unit vectors (after projections)...'
!    end if
!    !call clear_unitvecs(x)
!
!    if(mpi_cfg%myid==0) then
!      print*, 'Interpolation coords:  ',minval(self%zi),maxval(self%zi), &
!                                        minval(self%xi),maxval(self%xi), &
!                                        minval(self%yi),maxval(self%yi)
!      print*, 'Projection checking:  ',minval(self%proj_exp_e1),maxval(self%proj_exp_e1), &
!                                       minval(self%proj_exp_e2),maxval(self%proj_exp_e2), &
!                                       minval(self%proj_exp_e3),maxval(self%proj_exp_e3)
!    end if
!
!    self%flagcoordsi=.true.
!  end subroutine set_coordsi_neu3D


  !> do nothing stub - type extensions must override this to perform whatever load steps are needed for their data types
  subroutine load_size_neu3D(self)
    class(neutraldata3D_mpi), intent(inout) :: self

    return
  end subroutine load_size_neu3D


  !> do nothing stub
  subroutine load_grid_neu3D(self)
    class(neutraldata3D_mpi), intent(inout) :: self

    return
  end subroutine load_grid_neu3D


  subroutine load_data_neu3D(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutraldata3D_mpi), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
    integer :: lhorzn                        !number of horizontal grid points
    real(wp), dimension(:,:,:), allocatable :: paramall
    type(hdf5_file) :: hf
    character(:), allocatable :: fn

    UTsectmp = 0*t*dtmodel
    !! avoid unused argument warning

    lhorzn=self%lyn
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt,ymdtmp,UTsectmp)                !get the date for "next" params

    !read in the data from file
    if (mpi_cfg%myid==0) then    !root
      !in the 3D case we cannot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
      allocate(paramall(self%lzn,self%lxnall,self%lynall))     ! space to store a single neutral input parameter

      !print*, '  date and time (neutral3D):  ',ymdtmp,UTsectmp

      fn = date_filename(self%sourcedir,ymdtmp,UTsectmp) // ".h5"

      if (debug) print *, 'READ neutral 3D data from file: ',fn
      call hf%open(fn, action='r')

      call hf%read('/dn0all', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dnOall: non-finite value(s)'
      if (debug) print*, 'Min/max values for dnOall:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dnO,self%slabsizes,self%indx,self%dnO)
      call hf%read('/dnN2all', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dnN2all: non-finite value(s)'
      if (debug) print*, 'Min/max values for dnN2all:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dnN2,self%slabsizes,self%indx,self%dnN2)
      call hf%read('/dnO2all', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dnO2all: non-finite value(s)'
      if (debug) print*, 'Min/max values for dnO2all:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dnO2,self%slabsizes,self%indx,self%dnO2)
      call hf%read('/dTnall', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dTnall: non-finite value(s)'
      if (debug) print*, 'Min/max values for dTnall:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dTn,self%slabsizes,self%indx,self%dTn)
      call hf%read('/dvnrhoall', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnrhoall: non-finite value(s)'
      if (debug) print*, 'Min/max values for dvnrhoall:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dvnrho,self%slabsizes,self%indx,self%dvny)
      call hf%read('/dvnzall', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnzall: non-finite value(s)'
      if (debug) print*, 'Min/max values for dvnzall:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dvnz,self%slabsizes,self%indx,self%dvnz)
      call hf%read('/dvnxall', paramall)
      if (.not. all(ieee_is_finite(paramall))) error stop 'dvnxall: non-finite value(s)'
      if (debug) print*, 'Min/max values for dvnxall:  ',minval(paramall),maxval(paramall)
      call dneu_root2workers(paramall,tag%dvnx,self%slabsizes,self%indx,self%dvnx)

      call hf%close()
      deallocate(paramall)
    else     !workers
      !receive a subgrid copy of the data from root
      call dneu_workers_from_root(tag%dnO,self%dnO)
      call dneu_workers_from_root(tag%dnN2,self%dnN2)
      call dneu_workers_from_root(tag%dnO2,self%dnO2)
      call dneu_workers_from_root(tag%dTn,self%dTn)
      call dneu_workers_from_root(tag%dvnrho,self%dvny)
      call dneu_workers_from_root(tag%dvnz,self%dvnz)
      call dneu_workers_from_root(tag%dvnx,self%dvnx)
    end if


    if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
      print *, 'Min/max values for dnO:  ',mpi_cfg%myid,minval(self%dnO),maxval(self%dnO)
      print *, 'Min/max values for dnN:  ',mpi_cfg%myid,minval(self%dnN2),maxval(self%dnN2)
      print *, 'Min/max values for dnO2:  ',mpi_cfg%myid,minval(self%dnO2),maxval(self%dnO2)
      print *, 'Min/max values for dvnx:  ',mpi_cfg%myid,minval(self%dvnx),maxval(self%dvnx)
      print *, 'Min/max values for dvnrho:  ',mpi_cfg%myid,minval(self%dvny),maxval(self%dvny)
      print *, 'Min/max values for dvnz:  ',mpi_cfg%myid,minval(self%dvnz),maxval(self%dvnz)
      print *, 'Min/max values for dTn:  ',mpi_cfg%myid,minval(self%dTn),maxval(self%dTn)
    !  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
    end if
  end subroutine load_data_neu3D


  !> overriding procedure for updating neutral atmos (need additional rotation steps)
  subroutine update(self,cfg,dtmodel,t,x,ymd,UTsec)
    class(neutraldata3D_mpi), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel             ! need both model and input data time stepping
    real(wp), intent(in) :: t                   ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x            ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations

    ! execute a basic update
    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)

    ! FIXME: more efficient to rotate the winds only when interpolations are done...
    ! now we need to rotate velocity fields following interpolation (they are magnetic ENU prior to this step)
    call self%rotate_winds()

    ! FIXME: check if we need to further rotate these winds into geographic coordinates

    if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
      print*, ''
      print*, 'neutral data size:  ',mpi_cfg%myid,self%lzn,self%lxn,self%lyn
      print*, 'neutral data time:  ',ymd,UTsec
      print*, ''
      print *, 'Min/max values for dnOinext:  ',mpi_cfg%myid,minval(self%dnOinext),maxval(self%dnOinext)
      print *, 'Min/max values for dnNinext:  ',mpi_cfg%myid,minval(self%dnN2inext),maxval(self%dnN2inext)
      print *, 'Min/max values for dnO2inext:  ',mpi_cfg%myid,minval(self%dnO2inext),maxval(self%dnO2inext)
      print *, 'Min/max values for dvn1inext:  ',mpi_cfg%myid,minval(self%dvn1inext),maxval(self%dvn1inext)
      print *, 'Min/max values for dvn2inext:  ',mpi_cfg%myid,minval(self%dvn2inext),maxval(self%dvn2inext)
      print *, 'Min/max values for dvn3inext:  ',mpi_cfg%myid,minval(self%dvn3inext),maxval(self%dvn3inext)
      print *, 'Min/max values for dTninext:  ',mpi_cfg%myid,minval(self%dTninext),maxval(self%dTninext)
      print*, ''
      print *, 'Min/max values for dnOinow:  ',mpi_cfg%myid,minval(self%dnOinow),maxval(self%dnOinow)
      print *, 'Min/max values for dnNinow:  ',mpi_cfg%myid,minval(self%dnN2inow),maxval(self%dnN2inow)
      print *, 'Min/max values for dnO2inow:  ',mpi_cfg%myid,minval(self%dnO2inow),maxval(self%dnO2inow)
      print *, 'Min/max values for dvn1inow:  ',mpi_cfg%myid,minval(self%dvn1inow),maxval(self%dvn1inow)
      print *, 'Min/max values for dvn2inow:  ',mpi_cfg%myid,minval(self%dvn2inow),maxval(self%dvn2inow)
      print *, 'Min/max values for dvn3inow:  ',mpi_cfg%myid,minval(self%dvn3inow),maxval(self%dvn3inow)
      print *, 'Min/max values for dTninow:  ',mpi_cfg%myid,minval(self%dTninow),maxval(self%dTninow)
    end if
  end subroutine update


!  !> destructor for when object goes out of scope
!  subroutine destructor(self)
!    type(neutraldata3D) :: self
!
!    ! deallocate arrays from base inputdata class
!    call self%dissociate_pointers()
!
!    ! null pointers specific to parent neutraldata class
!    call self%dissociate_neutral_pointers()
!
!    ! now deallocate arrays specific to this extension
!    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
!    deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
!    deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
!    deallocate(self%ximat,self%yimat,self%zimat)
!
!    ! root has some extra data
!    if (mpi_cfg%myid==0) then
!      deallocate(self%extents,self%indx,self%slabsizes)
!      deallocate(self%xnall,self%ynall)
!    end if
!
!    ! set pointers to null
!    nullify(self%xi,self%yi,self%zi);
!    nullify(self%xn,self%yn,self%zn);
!    nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnx,self%dvny,self%dTn)
!  end subroutine destructor
end module neutraldata3Dobj_mpi
