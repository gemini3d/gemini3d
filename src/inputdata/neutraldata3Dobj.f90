module neutraldata3Dobj

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use meshobj, only: curvmesh
use config, only: gemini_cfg
use reader, only: get_simsize2,get_grid2,get_precip
use mpimod, only: mpi_integer,mpi_comm_world,mpi_status_ignore,mpi_realprec,mpi_cfg,tag=>gemini_mpi
use timeutils, only: dateinc,date_filename
use h5fortran, only: hdf5_file
use reader, only : get_simsize3
use pathlib, only: get_suffix,get_filename
use grid, only: gridflag

implicit none (type,external)
external :: mpi_send,mpi_recv
public :: neutraldata3D

!> type definition for 3D neutral data
type, extends(neutraldata) :: neutraldata3D
  ! source data coordinate pointers
  real(wp), dimension(:), pointer :: xn,yn,zn
  integer, pointer :: lxn,lyn,lzn
  real(wp), dimension(:), allocatable :: xnall,ynall
  integer :: lxnall,lynall

  ! work arrays needed by various procedures re: target coordinates
  real(wp), dimension(:,:,:), allocatable :: ximat,yimat,zimat
  real(wp), dimension(:), pointer :: zi,xi,yi

  ! source data pointers
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvnz,dvnx,dvny,dTn

  ! projection factors needed to rotate input data onto grid
  real(wp), dimension(:,:,:), allocatable :: proj_ezp_e1,proj_ezp_e2,proj_ezp_e3    
  real(wp), dimension(:,:,:), allocatable :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_exp_e1,proj_exp_e2,proj_exp_e3

  ! mpi-related information on subgrid extents and indices, only used on the root process; otherwise ignored
  real(wp), dimension(:,:), allocatable :: extents    ! min/max x,y,z of each worker
  integer, dimension(:,:), allocatable :: indx        ! indices for each workers' pieces of the neutral data
  integer, dimension(:,:), allocatable :: slabsizes  
  contains
    ! replacement for gridsize and gridload
    procedure :: load_sizeandgrid_neu3D
    procedure :: rotate_winds

    ! overriding procedures
    procedure :: update
    procedure :: init_storage

    ! bindings for deferred procedures
    procedure :: init=>init_neu3D
    procedure :: load_data=>load_data_neu3D
    procedure :: load_grid=>load_grid_neu3D    ! stub, does nothing see load_sizeandgrid_neu3D()
    procedure :: load_size=>load_size_neu3D    ! stub, does nothing "
    procedure :: set_coordsi=>set_coordsi_neu3D

    ! destructor
    final :: destructor
end type neutraldata3D


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
    class(neutraldata3D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel,dtdata
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation 
    integer :: lc1,lc2,lc3
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


  !> create storage for arrays needed specifically for 3D neutral input calculations, overrides the base class procedure
  subroutine init_storage(self)
    class(neutraldata3D), intent(inout) :: self
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

    ! input data coordinate arrays are set by load_gridandsize()

    ! allocate target coords, for neutral3D the standard set (coord1i, etc.) are done in set_coordsi()
    allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
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
  end subroutine init_storage


  !> do nothing stub
  subroutine load_size_neu3D(self)
    class(neutraldata3D), intent(inout) :: self

  end subroutine load_size_neu3D


  !> do nothing stub
  subroutine load_grid_neu3D(self)
    class(neutraldata3D), intent(inout) :: self

  end subroutine load_grid_neu3D


  !> load source data size and grid information and communicate to worker processes.  
  !    Note that this routine will allocate sizes for source coordinates grids in constrast 
  !    with other inputdata type extensions which have separate load_size, allocate, and 
  !    load_grid procedures.  
  subroutine load_sizeandgrid_neu3D(self,cfg)
    class(neutraldata3D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), dimension(:), allocatable :: xn,yn             ! for root to break off pieces of the entire grid array
    integer :: ix1,ix2,ix3,ihorzn,izn,iid,ierr
    integer :: lxntmp,lyntmp                                   ! local copies for root, eventually these need to be stored in object
    real(wp) :: maxzn
    real(wp), dimension(2) :: xnrange,ynrange                ! these eventually get stored in extents
    integer, dimension(6) :: indices                         ! these eventually get stored in indx
    integer :: ixn,iyn
    integer :: lxn,lyn
    real(wp) :: meanxn,meanyn

    if (mpi_cfg%myid==0) then    !root must establish the size of the grid based on input file and distribute to workers
      print '(A,/,A)', 'READ neutral size from:', self%sourcedir
      call get_simsize3(self%sourcedir, lx1=self%lxnall, lx2all=self%lynall, lx3all=self%lzn)
      print *, 'Neutral data has lx,ly,lz size:  ',self%lxnall,self%lynall,self%lzn, &
                   ' with spacing dx,dy,dz',cfg%dxn,cfg%drhon,cfg%dzn
      if (self%lxnall < 1 .or. self%lynall < 1 .or. self%lzn < 1) then
        write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
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
        call mpi_send(self%lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
        call mpi_send(self%zn,self%lzn,mpi_realprec,iid,tag%zn,MPI_COMM_WORLD,ierr)
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
        call mpi_recv(xnrange,2,mpi_realprec,iid,tag%xnrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        call mpi_recv(ynrange,2,mpi_realprec,iid,tag%ynrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
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
        call mpi_send(lyn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
        call mpi_send(lxn,1,MPI_INTEGER,iid,tag%lx,MPI_COMM_WORLD,ierr)
        allocate(xn(lxn),yn(lyn))
        xn=self%xnall(self%indx(iid,3):self%indx(iid,4))
        yn=self%ynall(self%indx(iid,5):self%indx(iid,6))
        call mpi_send(xn,lxn,mpi_realprec,iid,tag%xn,MPI_COMM_WORLD,ierr)
        call mpi_send(yn,lyn,mpi_realprec,iid,tag%yn,MPI_COMM_WORLD,ierr)
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
      call mpi_recv(self%lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

      ! allocate space for target coordinate and bind alias
      allocate(self%coord1(self%lzn))
      self%zn=>self%coord1

      ! receive data from root
      call mpi_recv(self%zn,self%lzn,mpi_realprec,0,tag%zn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      maxzn=maxval(self%zn)
    
      !calculate the extent of my grid
      call slabrange(maxzn,self%ximat,self%yimat,self%zimat,cfg%sourcemlat,xnrange,ynrange,gridflag)
    
      !send ranges to root
      call mpi_send(xnrange,2,mpi_realprec,0,tag%xnrange,MPI_COMM_WORLD,ierr)
      call mpi_send(ynrange,2,mpi_realprec,0,tag%ynrange,MPI_COMM_WORLD,ierr)
    
      !receive my sizes from root, allocate then receive my pieces of the grid
      call mpi_recv(self%lxn,1,MPI_INTEGER,0,tag%lx,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call mpi_recv(self%lyn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

      ! at this point we can allocate space for the source coordinates and bind aliases as needed
      allocate(self%coord2(self%lxn),self%coord3(self%lyn))
      self%xn=>self%coord2; self%yn=>self%coord3;        ! input data coordinates

      ! recieve data from root
      call mpi_recv(self%xn,self%lxn,mpi_realprec,0,tag%xn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call mpi_recv(self%yn,self%lyn,mpi_realprec,0,tag%yn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end if

    self%flagdatasize=.true.
  end subroutine load_sizeandgrid_neu3D


  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
  subroutine set_coordsi_neu3D(self,cfg,x)
    class(neutraldata3D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
    real(wp) :: xp,yp
    real(wp), dimension(3) :: ezp,eyp,tmpvec,exprm
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid,ierr


    ! Space for coordinate sites and projections in neutraldata3D object
    allocate(self%coord1i(x%lx1*x%lx2*x%lx3),self%coord2i(x%lx1*x%lx2*x%lx3),self%coord3i(x%lx1*x%lx2*x%lx3))
    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i;     ! coordinates of interpolation sites
    allocate(self%ximat(x%lx1,x%lx2,x%lx3),self%yimat(x%lx1,x%lx2,x%lx3),self%zimat(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_ezp_e1(x%lx1,x%lx2,x%lx3),self%proj_ezp_e2(x%lx1,x%lx2,x%lx3),self%proj_ezp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_eyp_e1(x%lx1,x%lx2,x%lx3),self%proj_eyp_e2(x%lx1,x%lx2,x%lx3),self%proj_eyp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_exp_e1(x%lx1,x%lx2,x%lx3),self%proj_exp_e2(x%lx1,x%lx2,x%lx3),self%proj_exp_e3(x%lx1,x%lx2,x%lx3)) 

    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2 - cfg%sourcemlat*pi/180
    
    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
    if (mpi_cfg%myid==0) then
      print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
    end if

    self%zimat=x%alt     !vertical coordinate is just altitude array already stored in grid object
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          ! interpolation based on geomag
          theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle

          !print*, ' center NS set',shape(self%zi),shape(self%zimat),theta2,x%theta(ix1,ix2,ix3)

          if (x%lx2/=1) then
            phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
          else
            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
          end if
 
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
    
          ! coordinates from distances
          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
            yp= -yp            !do we want an abs here to be safe
          end if
          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
            xp= -xp
          end if

          self%ximat(ix1,ix2,ix3)=xp     !eastward distance
          self%yimat(ix1,ix2,ix3)=yp     !northward distance
    
          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
          ezp=x%er(ix1,ix2,ix3,:)
    
          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca
    
          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca
    
          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca
    
          eyp= -x%etheta(ix1,ix2,ix3,:)
    
          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca
    
          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca
    
          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca
    
          exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)
    
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
    
    !Assign values for flat lists of grid points
    if (mpi_cfg%myid==0) then
      print*, '...Packing interpolation target points...'
    end if
    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
    self%yi=pack(self%yimat,.true.)
    self%xi=pack(self%ximat,.true.)

    ! FIXME: do we need to have the new grid code clear its unit vectors?  Or maybe this isn't a huge waste of memory???
    if (mpi_cfg%myid==0) then
      print*, '...Clearing out unit vectors (after projections)...'
    end if
    !call clear_unitvecs(x)
    
    if(mpi_cfg%myid==0) then
      print*, 'Interpolation coords:  ',minval(self%zi),maxval(self%zi), &
                                        minval(self%xi),maxval(self%xi), &
                                        minval(self%yi),maxval(self%yi)
      print*, 'Projection checking:  ',minval(self%proj_exp_e1),maxval(self%proj_exp_e1), &
                                       minval(self%proj_exp_e2),maxval(self%proj_exp_e2), &
                                       minval(self%proj_exp_e3),maxval(self%proj_exp_e3)
    end if
    
    self%flagcoordsi=.true.
  end subroutine set_coordsi_neu3D


  subroutine load_data_neu3D(self,t,dtmodel,ymdtmp,UTsectmp)
    class(neutraldata3D), intent(inout) :: self
    real(wp), intent(in) :: t,dtmodel
    integer, dimension(3), intent(inout) :: ymdtmp
    real(wp), intent(inout) :: UTsectmp
    integer :: iid,ierr
    integer :: lhorzn                        !number of horizontal grid points
    real(wp), dimension(:,:,:), allocatable :: paramall
    type(hdf5_file) :: hf
    character(:), allocatable :: fn
        
    lhorzn=self%lyn
    ymdtmp = self%ymdref(:,2)
    UTsectmp = self%UTsecref(2)
    call dateinc(self%dt,ymdtmp,UTsectmp)                !get the date for "next" params

    !read in the data from file
    if (mpi_cfg%myid==0) then    !root
      !in the 3D case we cannot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
      allocate(paramall(self%lzn,self%lxnall,self%lynall))     ! space to store a single neutral input parameter

      print*, '  date and time (neutral3D):  ',ymdtmp,UTsectmp

      fn=date_filename(self%sourcedir,ymdtmp,UTsectmp)
      fn=get_filename(fn)
      if (debug) print *, 'READ neutral 3D data from file: ',fn
      if (get_suffix(fn)=='.h5') then
        call hf%open(fn, action='r')
      else
        error stop '3D neutral input only supported for hdf5 files; please regenerate input'
      end if
    
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
    class(neutraldata3D), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: dtmodel             ! need both model and input data time stepping
    real(wp), intent(in) :: t                   ! simulation absoluate time for which perturabation is to be computed
    class(curvmesh), intent(in) :: x            ! mesh object
    integer, dimension(3), intent(in) :: ymd    ! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec               ! UT seconds for which we with to compute perturbations

    ! execute a basic update
    call self%update_simple(cfg,dtmodel,t,x,ymd,UTsec)

    ! now we need to rotate velocity fields following interpolation (they are magnetic ENU prior to this step)
    call self%rotate_winds()

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


  !> This subroutine takes winds stored in self%dvn?inow and applies a rotational transformation onto the 
  !      grid object for this simulation
  subroutine rotate_winds(self)
    class(neutraldata3D), intent(inout) :: self
    integer :: ix1,ix2,ix3
    real(wp) :: vnx,vny,vnz

    ! do rotations one grid point at a time to cut down on temp storage needed
    do ix3=1,self%lc3i
      do ix2=1,self%lc2i
        do ix1=1,self%lc3i
          vnz=self%dvn1inext(ix1,ix2,ix3)
          vnx=self%dvn2inext(ix1,ix2,ix3)
          vny=self%dvn3inext(ix1,ix2,ix3)
          self%dvn1inext(ix1,ix2,ix3)=vnz*self%proj_ezp_e1(ix1,ix2,ix3) + vnx*self%proj_exp_e1(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e1(ix1,ix2,ix3)
          self%dvn2inext(ix1,ix2,ix3)=vnz*self%proj_ezp_e2(ix1,ix2,ix3) + vnx*self%proj_exp_e2(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e2(ix1,ix2,ix3)
          self%dvn3inext(ix1,ix2,ix3)=vnz*self%proj_ezp_e3(ix1,ix2,ix3) + vnx*self%proj_exp_e3(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e3(ix1,ix2,ix3)
        end do
      end do
    end do
  end subroutine rotate_winds


  !> destructor for when object goes out of scope
  subroutine destructor(self)
    type(neutraldata3D) :: self

    ! deallocate arrays from base inputdata class
    call self%dissociate_pointers()

    ! null pointers specific to parent neutraldata class
    call self%dissociate_neutral_pointers()

    ! now deallocate arrays specific to this extension
    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
    deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
    deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
    deallocate(self%extents,self%indx,self%slabsizes)
    deallocate(self%ximat,self%yimat,self%zimat)

    ! root has some extra data
    if (mpi_cfg%myid==0) then
      deallocate(self%xnall,self%ynall)
    end if

    ! set pointers to null
    nullify(self%xi,self%yi,self%zi);
    nullify(self%xn,self%yn,self%zn);
    nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnx,self%dvny,self%dTn)
  end subroutine destructor
end module neutraldata3Dobj
