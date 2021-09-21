module neutraldata2Daxisymmobj

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use meshobj, only: curvmesh
use config, only: gemini_cfg
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use neutraldata2Dobj, only: neutraldata2D
use reader, only: get_simsize3
use mpimod, only: mpi_integer,mpi_comm_world,mpi_status_ignore,mpi_realprec,mpi_cfg,tag=>gemini_mpi

implicit none (type, external)
external :: mpi_send,mpi_recv
public :: neutraldata2Daxisymm

!> type extension for neutral 2D axisymmetric input data
type, extends(neutraldata2D) :: neutraldata2Daxisymm
  integer, pointer :: lrhon
  real(wp), dimension(:), pointer :: rhon
  real(wp), dimension(:), pointer :: rhoi

  contains
    procedure :: init=>init_neu2Daxisymm
    procedure :: load_sizeandgrid_neu2D=>load_sizeandgrid_neu2Daxisymm
    procedure :: set_coordsi=>set_coordsi_neu2Daxisymm
    final :: destructor
end type neutraldata2Daxisymm

contains
  subroutine init_neu2Daxisymm(self,cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)
    class(neutraldata2Daxisymm), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    character(*), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x
    real(wp), intent(in) :: dtmodel,dtdata
    integer, dimension(3), intent(in) :: ymd            ! target date of initiation
    real(wp), intent(in) :: UTsec                       ! target time of initiation 
    character(:), allocatable :: strname

    ! basic init for any 2D neutral input
    call self%init_neu2D_simple(cfg,sourcedir,x,dtmodel,dtdata,ymd,UTsec)

    ! append type of interp. to dataname
    strname=self%dataname//' axisymmetric'     ! append type of 2D interpolation to name
    call self%set_name(strname)                ! overwrite generic neutral 2D data name

    ! bind axisymmetric specific pointers for convenience, in case they are needed elsewhere
    self%lrhon=>self%lhorzn
    self%rhon=>self%horzn
    self%rhoi=>self%horzi
  end subroutine init_neu2Daxisymm


  !! FIXME:  currently hardcoded for axisymmetric coords.  Needs to be specific to coordinate system.  
  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
  subroutine set_coordsi_neu2Daxisymm(self,cfg,x)
    class(neutraldata2Daxisymm), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
    real(wp) :: xp,yp
    real(wp), dimension(3) :: ezp,erhop,tmpvec,exprm
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid,ierr


    ! Space for coordinate sites and projections in neutraldata2D object
    allocate(self%coord1i(x%lx1*x%lx2*x%lx3),self%coord3i(x%lx1*x%lx2*x%lx3))
    self%zi=>self%coord1i; self%horzi=>self%coord3i;     ! coordinates of interpolation sites
    allocate(self%horzimat(x%lx1,x%lx2,x%lx3),self%zimat(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_ezp_e1(x%lx1,x%lx2,x%lx3),self%proj_ezp_e2(x%lx1,x%lx2,x%lx3),self%proj_ezp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_ehorzp_e1(x%lx1,x%lx2,x%lx3),self%proj_ehorzp_e2(x%lx1,x%lx2,x%lx3),self%proj_ehorzp_e3(x%lx1,x%lx2,x%lx3))

    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2-cfg%sourcemlat*pi/180
    
    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
    if (mpi_cfg%myid==0) then
      print *, 'Computing alt,radial distance values for plasma grid and completing rotations'
    end if
    self%zimat=x%alt     !vertical coordinate
    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          !INTERPOLATION BASED ON GEOMAGNETIC COORDINATES
          theta2=x%theta(ix1,ix2,ix3)                    !field point zenith angle
          if (x%lx2/=1 .and. x%lx3/=1) then
            phi2=x%phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
          else
            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
          end if
    
          !COMPUTE DISTANCES
          gammarads=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)     !this is actually cos(gamma)
          if (gammarads > 1) then     !handles weird precision issues in 2D
            gammarads = 1
          else if (gammarads < -1) then
            gammarads= -1
          end if
          gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
          self%horzimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''
    
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
!            gamma2 = 1
!          else if (gamma2< -1) then
!            gamma2= -1
!          end if
!          gamma2=acos(gamma2)
!    
!          xp=Re*gamma1
!          yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement
!    
!          !COMPUTE COORDINATES FROM DISTANCES
!          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
!            yp = -yp            !do we want an abs here to be safe
!          end if
!          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
!            xp = -xp
!          end if
!          phip=atan2(yp,xp)
!    
!          if(flagcart) then
!            yimat(ix1,ix2,ix3)=yp
!          end if
    
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
    
          !if (flagcart) then
          !  eyp= -x%etheta(ix1,ix2,ix3,:)
          !
          !  tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          !  tmpsca=sum(tmpvec)
          !  proj_eyp_e1(ix1,ix2,ix3)=tmpsca
          !
          !  tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          !  tmpsca=sum(tmpvec)
          !  proj_eyp_e2(ix1,ix2,ix3)=tmpsca
          !
          !  tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          !  tmpsca=sum(tmpvec)
          !  proj_eyp_e3(ix1,ix2,ix3)=tmpsca
          !else
            erhop=cos(phip)*x%e3(ix1,ix2,ix3,:) - sin(phip)*x%etheta(ix1,ix2,ix3,:)     !unit vector for azimuth (referenced from epicenter - not geocenter!!!) in cartesian geocentric-geomagnetic coords.
    
            tmpvec=erhop*x%e1(ix1,ix2,ix3,:)
            tmpsca=sum(tmpvec)
            self%proj_ehorzp_e1(ix1,ix2,ix3)=tmpsca
    
            tmpvec=erhop*x%e2(ix1,ix2,ix3,:)
            tmpsca=sum(tmpvec)
            self%proj_ehorzp_e2(ix1,ix2,ix3)=tmpsca
    
            tmpvec=erhop*x%e3(ix1,ix2,ix3,:)
            tmpsca=sum(tmpvec)
            self%proj_ehorzp_e3(ix1,ix2,ix3)=tmpsca
          !end if
        end do
      end do
    end do
    
    !Assign values for flat lists of grid points
    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation ffunctions
    !if (flagcart) then
    !  yi=pack(yimat,.true.)
    !else
      self%horzi=pack(self%horzimat,.true.)
    !end if
    
    !GRID UNIT VECTORS NO LONGER NEEDED ONCE PROJECTIONS ARE CALCULATED...
    !call clear_unitvecs(x)
    
    !PRINT OUT SOME BASIC INFO ABOUT THE GRID THAT WE'VE LOADED
    if (mpi_cfg%myid==0 .and. debug) then
      !if (flagcart) then
      !  print *, 'Min/max yn,zn values',minval(yn),maxval(yn),minval(zn),maxval(zn)
      !  print *, 'Min/max yi,zi values',minval(yi),maxval(yi),minval(zi),maxval(zi)
      !else
        print *, 'Min/max rhon,zn values',minval(self%horzn),maxval(self%horzn),minval(self%zn),maxval(self%zn)
        print *, 'Min/max rhoi,zi values',minval(self%horzi),maxval(self%horzi),minval(self%zi),maxval(self%zi)
      !end if
    
      print *, 'Source lat/long:  ',cfg%sourcemlat,cfg%sourcemlon
      print *, 'Plasma grid lat range:  ',minval(x%glat(:,:,:)),maxval(x%glat(:,:,:))
      print *, 'Plasma grid lon range:  ',minval(x%glon(:,:,:)),maxval(x%glon(:,:,:))
    end if

    self%flagcoordsi=.true.
  end subroutine set_coordsi_neu2Daxisymm


  !! FIXME: may be specific to axisymmetric vs. cartesian
  !> load source data size and grid information and communicate to worker processes.  
  !    Note that this routine will allocate sizes for source coordinates grids in constrast 
  !    with other inputdata type extensions which have separate load_size, allocate, and 
  !    load_grid procedures.  
  subroutine load_sizeandgrid_neu2Daxisymm(self,cfg)
    class(neutraldata2Daxisymm), intent(inout) :: self
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
    real(wp) :: dhorzn

    !horizontal grid spacing
    dhorzn=cfg%drhon
    
    !Establish the size of the grid based on input file and distribute to workers
    if (mpi_cfg%myid==0) then    !root
      print '(A,/,A)', 'Inputting neutral size from:  ',self%sourcedir
    
    ! bit of a tricky issue here; for neutral input, according to makedneuframes.m, the first integer in the size file is
    !  the horizontal grid point count for the input - which get_simsize3 interprets as lx1...
    call get_simsize3(cfg%sourcedir, lx1=self%lhorzn, lx2all=self%lzn)
    self%lxn=1     ! treat as a 3D dataset with singleton dimension along x
    
      print *, 'Neutral data has lhorzn,lz size:  ',self%lhorzn,self%lzn,' with spacing dhorzn,dz',dhorzn,cfg%dzn
      if (self%lhorzn < 1 .or. self%lzn < 1) then
        write(stderr,*) 'ERROR: reading ' // self%sourcedir
        error stop 'neutral:gridproj_dneu2D: grid size must be strictly positive'
      endif
      do iid=1,mpi_cfg%lid-1
        call mpi_send(self%lhorzn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
        call mpi_send(self%lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
      end do
    else                 !workers
      call mpi_recv(self%lhorzn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call mpi_recv(self%lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end if
    
    !Everyone must allocate space for the grid of input data
    allocate(self%zn(self%lzn))    !these are module-scope variables
    !if (flagcart) then
    !  allocate(rhon(1))  !not used in Cartesian code so just set to something
    allocate(self%horzn(self%lhorzn))    ! FIXME: default to axisymmetric?
    !  lyn=lhorzn
    !else
    !  allocate(rhon(lhorzn))
    !  allocate(yn(1))    !not used in the axisymmetric code so just initialize to something
    !  lrhon=lhorzn
    !end if
    
    !Note that the second dimension ("longitude") is singleton so that we are able to also use these vars for 3D input
    !allocate(dnO(lzn,1,lhorzn),dnN2(lzn,1,lhorzn),dnO2(lzn,1,lhorzn),dvnrho(lzn,1,lhorzn),dvnz(lzn,1,lhorzn),dTn(lzn,1,lhorzn))
    
    !Define a grid (input data) by assuming that the spacing is constant
    !if (flagcart) then     !Cartesian neutral simulation
    !  yn=[ ((real(ihorzn, wp)-1)*dhorzn, ihorzn=1,lhorzn) ]
    !  meanyn=sum(yn,1)/size(yn,1)
    !  yn=yn-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
    !else
    self%horzn=[ ((real(ihorzn, wp)-1)*dhorzn, ihorzn=1,self%lhorzn) ]
    !end if
    self%zn=[ ((real(izn, wp)-1)*cfg%dzn, izn=1,self%lzn) ]
    
    if (mpi_cfg%myid==0) then
    !  if (flagcart) then
    !    print *, 'Creating neutral grid with y,z extent:',minval(yn),maxval(yn),minval(zn),maxval(zn)
    !  else
        print *, 'Creating neutral grid with rho,z extent:  ',minval(self%horzn),maxval(self%horzn),minval(self%zn),maxval(self%zn)
    !  end if
    end if

    self%flagdatasize=.true.
  end subroutine load_sizeandgrid_neu2Daxisymm

  subroutine destructor(self)
    type(neutraldata2Daxisymm), intent(inout) :: self

    ! de facto neutral2D destructor takes care of most everything
    call self%dissociate_neutral2D_pointers()

    ! now extension-specific quantities
    nullify(self%lrhon,self%rhon,self%rhoi)
  end subroutine destructor
end module neutraldata2Daxisymmobj
