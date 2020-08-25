module neutral

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit

use reader, only : get_simsize3, get_neutral2, get_neutral3
use phys_consts, only: wp, lnchem, pi, re, debug
use grid, only: lx1, lx2, lx3, clear_unitvecs, gridflag
use mesh, only: curvmesh
use interpolation, only : interp2, interp3
use timeutils, only : doy_calc,dateinc, date_filename, find_lastdate
use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
myid, lid, mpi_realprec, tag=>gemini_mpi
use config, only: gemini_cfg

! also links gtd7 from vendor/msis00/

implicit none (type, external)
private
public :: Tnmsis, neutral_atmos, make_dneu, clear_dneu, neutral_perturb, neutral_update, init_neutrals


interface ! atmos.f90
module subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,nn,Tn,vn1,vn2,vn3)
integer, intent(in) :: ymd(3)
real(wp), intent(in) :: UTsecd
real(wp), dimension(:,:,:), intent(in) :: glat,glon,alt
real(wp), intent(in) :: activ(3)
real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3),lnchem), intent(out) :: nn
real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(out) :: Tn
real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)), intent(out) :: vn1,vn2,vn3
end subroutine neutral_atmos
end interface

external :: mpi_send, mpi_recv

!! ALL ARRAYS THAT FOLLOW ARE USED WHEN INCLUDING NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
!! ARRAYS TO STORE THE NEUTRAL GRID INFORMATION
!! as long as the neutral module is in scope these persist and do not require a "save"; this variable only used by the axisymmetric interpolation
real(wp), dimension(:), allocatable, private :: rhon     !used for axisymmetric 2D simulations
real(wp), dimension(:), allocatable, private :: yn    !used in cartesian 2D and 3D interpolation
real(wp), dimension(:), allocatable, private :: zn
real(wp), dimension(:), allocatable, private :: xn    !for 3D cartesian interpolation
integer, private :: lrhon,lzn,lyn,lxn


!! STORAGE FOR NEUTRAL SIMULATION DATA.
! These will be singleton in the second dimension (longitude) in the case of 2D interpolation...
!! THESE ARE INCLUDED AS MODULE VARIATIONS TO AVOID HAVING TO REALLOCATE AND DEALLOCIATE EACH TIME WE NEED TO INTERP
real(wp), dimension(:,:,:), allocatable, private :: dnO,dnN2,dnO2,dvnrho,dvnz,dvnx,dTn


!!full grid parameters for root to store input from files.
real(wp), dimension(:), allocatable, private :: xnall
real(wp), dimension(:), allocatable, private :: ynall
integer, private :: lxnall,lynall

real(wp), dimension(:,:,:), allocatable, private :: dnOall,dnN2all,dnO2all,dvnrhoall,dvnzall,dvnxall,dTnall


!ARRAYS TO STORE NEUTRAL DATA THAT HAS BEEN INTERPOLATED
real(wp), dimension(:,:,:), allocatable, private :: dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev, &
                                                   dvn1iprev,dvn2iprev,dvn3iprev,dvnxiprev
real(wp), private :: tprev
integer, dimension(3), private :: ymdprev
!! time corresponding to "prev" interpolated data

real(wp), private :: UTsecprev
real(wp), dimension(:,:,:), allocatable, private :: dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext, &
                                                   dTninext,dvn1inext,dvn2inext,dvn3inext,dvnxinext
real(wp), private :: tnext
integer, dimension(3), private :: ymdnext
real(wp), private :: UTsecnext

!! data at current time level, (centered in time between current time step and next)
real(wp), dimension(:,:,:), allocatable, protected :: dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow


!SPACE TO STORE PROJECTION FACTORS
real(wp), dimension(:,:,:), allocatable, private :: proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2,proj_erhop_e3,proj_ezp_e3    !these projections are used in the axisymmetric interpolation
real(wp), dimension(:,:,:), allocatable, private :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3    !these are for Cartesian projections
real(wp), dimension(:,:,:), allocatable, private :: proj_exp_e1,proj_exp_e2,proj_exp_e3


!PLASMA GRID ZI AND RHOI LOCATIONS FOR INTERPOLATIONS
real(wp), dimension(:), allocatable, private :: zi,yi,xi,rhoi    !this is to be a flat listing of sites on the, rhoi only used in axisymmetric and yi only in cartesian


!USED FOR 3D INTERPOLATION WHERE WORKER DIVISIONS ARE COMPLICATED (note that the first dim starts at zero so it matches mpi ID)
real(wp), dimension(:,:), private, allocatable :: extents    !roots array that is used to store min/max x,y,z of each works
integer, dimension(:,:), private, allocatable :: indx       !roots array that contain indices for each workers needed piece of the neutral data
integer, dimension(:,:), private, allocatable :: slabsizes


!! BASE MSIS ATMOSPHERIC STATE ON WHICH TO APPLY PERTURBATIONS
real(wp), dimension(:,:,:,:), allocatable, protected :: nnmsis
real(wp), dimension(:,:,:), allocatable, protected :: Tnmsis
real(wp), dimension(:,:,:), allocatable, protected :: vn1base,vn2base,vn3base

contains


subroutine init_neutrals(dt,t,cfg,ymd,UTsec,x,nn,Tn,vn1,vn2,vn3)

!> initializes neutral atmosphere by:
!    1)  allocating storage space
!    2)  establishing initial background
!    3)  priming file input so that we have an initial perturbed state to start from (necessary for restart)

real(wp), intent(in) :: dt,t
type(gemini_cfg), intent(in) :: cfg
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
type(curvmesh), intent(inout) :: x    ! unit vecs may be deallocated after first setup
real(wp), dimension(:,:,:,:), intent(out) :: nn
real(wp), dimension(:,:,:), intent(out) :: Tn
real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

real(wp) :: tstart,tfin


!! allocation neutral module scope variables so there is space to store all the file input and do interpolations
call make_dneu()

!! call msis to get an initial neutral background atmosphere
if (myid==0) call cpu_time(tstart)
call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,nn,Tn,vn1,vn2,vn3)
if (myid==0) then
  call cpu_time(tfin)
  print *, 'Initial neutral background at time:  ',ymd,UTsec,' calculated in time:  ',tfin-tstart
end if

if (cfg%flagdneu==1) then
  !! find the last input data preceding the milestone/initial condition that we start with
  call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,cfg%dtneu,ymdtmp,UTsectmp)

  !! Loads the neutral input file corresponding to the "first" time step of the simulation to prevent the first interpolant
  !  from being zero and causing issues with restart simulations.  I.e. make sure the neutral buffers are primed for restart
  !  This requires us to load file input twice, once corresponding to the initial frame and once for the "first, next" frame.
  tprev=UTsectmp-UTsec-2._wp*cfg%dtneu
  tnext=tprev+cfg%dtneu
  if (myid==0) print*, '!!!Attempting initial load of neutral dynamics files!!!' // &
                           ' This is a workaround to insure compatibility with restarts...',ymdtmp,UTsectmp
  !! We essentially are loading up the data corresponding to halfway betwween -dtneu and t0 (zero).  This will load
  !   two time levels back so when tprev is incremented twice it will be the true tprev corresponding to first time step
  call neutral_perturb(cfg,dt,cfg%dtneu,tnext+cfg%dtneu/2._wp,ymdtmp,UTsectmp-cfg%dtneu, &
                        x,nn,Tn,vn1,vn2,vn3)  !abs time arg to be < 0

  if (myid==0) print*, 'Now loading initial next file for neutral perturbations...'
  !! Now compute perturbations for the present time (zero), this moves the primed variables in next into prev and then
  !  loads up a current state so that we get a proper interpolation for the first time step.
  call neutral_perturb(cfg,dt,cfg%dtneu,0._wp,ymdtmp,UTsectmp,x,nn,Tn,vn1,vn2,vn3)    !t-dt so we land exactly on start time
end if
end subroutine init_neutrals


subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,nn,Tn,vn1,vn2,vn3)

!THIS IS  WRAPPER FOR THE NEUTRAL PERTURBATION CODES THAT DO EITHER
!AXISYMMETRIC OR CARTESIAN OR 3D INTERPOLATION

type(gemini_cfg), intent(in) :: cfg
real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd     !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec

type(curvmesh), intent(inout) :: x                !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3


if (cfg%interptype==0) then                           !cartesian interpolation drho inputs (radial distance) will be interpreted as dy (horizontal distance)
!  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,neudir,drhon,dzn,meanlat,meanlong,x,nn,Tn,vn1,vn2,vn3)
  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else if (cfg%interptype==1) then                      !axisymmetric interpolation
  call neutral_perturb_axisymm(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else if (cfg%interptype==3) then                      !3D interpolation drhon is takent to be dyn (northward distance)
  call neutral_perturb_3D(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else
  error stop '...Invalid interpolation type specified from input file...'
end if

end subroutine neutral_perturb


!FOR CONSISTENCY I'D LIKE TO STRUCTURE NEUTRAL PERTURB OPERATIONS LIKE GRAVITY IS HANDLED IN THE GRID MODULE, I.E. HAVE AN EXPLICIT CONSTRUCTORS/DESTRUCTOR TYPE ROUTINE THAT HANDLES ALLOCATION AND DEALLOCATION, WHICH WILL CLEAN UP THE NEUTRAL_PERTURB SUBROUTINE, I.E. REMOVE ALLOCATES OF PERSISTENT MODULE VARIABLES.
subroutine neutral_perturb_axisymm(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CYLINDRICAL COORDINATES.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid!,irhon,izn
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2._wp>=tnext .or. t<0._wp) then   !negative time means that we need to load the first frame

  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load??? or just assume everything starts at zero?  This needs to somehow check for an existing file under certain conditiosn, maybe if it==1???  Actually we don't even need that we can just check that the neutral grid is allocated (or not)
    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    call gridproj_dneu2D(cfg,.false.,x)    !set false to denote not Cartesian; here and below...
  end if

  !Read in neutral data from a file
  call read_dneu2D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp,.false.)

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
  end if
  call spaceinterp_dneu2D(.false.)

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if !done loading frame data...

!Interpolation in time
call timeinterp_dneu(t,dt,dnOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!Add interpolated perturbations to reference atmosphere arrays
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_axisymm


!! THIS SHARES SO MUCH CODE WITH THE AXISYMMETRIC VERSION THAT THEY SHOULD PROBABLY BE COMBINED
subroutine neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER.
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CARTESIAN COORDINATES.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2d0 >= tnext .or. t <= 0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    tprev=t
    tnext=t+cfg%dtneu

    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    call gridproj_dneu2D(cfg,.true.,x)    !set true to denote Cartesian...
  end if

  !Read in neutral data from a file
  call read_dneu2D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp,.true.)

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
  end if
  call spaceinterp_dneu2D(.true.)

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if

!Interpolation in time
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!> update neutral atmosphere
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_cart


subroutine neutral_perturb_3D(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER.
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CARTESIAN COORDINATES AND IN THREE DIMENSIONS.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

real(wp) :: starttime,endtime


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2d0>=tnext .or. t<=0d0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    tprev=t
    tnext=t+cfg%dtneu

    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    if (myid==0 .and. debug) then
      print*, 'Creating a neutral grid...'
    end if
    call gridproj_dneu3D(cfg,x)
  end if

  !Read in neutral data from a file
  if (myid==0 .and. debug) then
    print*, 'Reading in data from neutral file'
    call cpu_time(starttime)
  end if
  call read_dneu3D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp)
  if (myid==0 .and. debug) then
    call cpu_time(endtime)
    print*, 'Neutral data input required time:  ',endtime-starttime
  end if

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
    call cpu_time(starttime)
  end if
  call spaceinterp_dneu3D()
  if (myid==0 .and. debug) then
    call cpu_time(endtime)
    print*, 'Spatial interpolation in 3D took time:  ',endtime-starttime
  end if

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if

!Interpolation in time
if (myid==0 .and. debug) then
  print*, 'Interpolating in time'
end if
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!> update neutral atmosphere
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_3D


subroutine gridproj_dneu2D(cfg,flagcart,x)

!Read in the grid for the neutral data and project unit vectors into the appropriiate directions.
!Also allocate module-scope variables for storing neutral perturbations read in from input files.

type(gemini_cfg), intent(in) :: cfg
logical, intent(in) :: flagcart              !whether or not the input data are to be interpreted as Cartesian
type(curvmesh), intent(inout) :: x           !inout to allow deallocation of unit vectors once we are done with them, should consider exporting this to another functino to be called from main program to avoid having x writeable...

real(wp) :: dhorzn           !neutral grid spacing in horizontal "rho or y" and vertical directions
integer :: lhorzn
real(wp) :: meanyn

real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
real(wp) :: xp,yp
real(wp), dimension(3) :: erhop,ezp,eyp,tmpvec
real(wp) :: tmpsca

integer :: ix1,ix2,ix3,ihorzn,izn,iid,ierr
real(wp), dimension(x%lx1,x%lx2,x%lx3) :: zimat,rhoimat,yimat


!horizontal grid spacing
dhorzn=cfg%drhon

!Establish the size of the grid based on input file and distribute to workers
if (myid==0) then    !root
  print '(A,/,A)', 'Inputting neutral size from:  ',cfg%sourcedir

  ! bit of a tricky issue here; for neutral input, according to makedneuframes.m, the first integer in the size file is
  !  the horizontal grid point count for the input - which get_simsize3 interprets as lx1...
  !call get_simsize3(cfg%sourcedir, lx1=lzn, lx2all=lhorzn)
call get_simsize3(cfg%sourcedir, lx1=lhorzn, lx2all=lzn)

  print *, 'Neutral data has lhorzn,lz size:  ',lhorzn,lzn,' with spacing dhorzn,dz',dhorzn,cfg%dzn
  if (lhorzn < 1 .or. lzn < 1) then
    write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
    error stop 'neutral:gridproj_dneu2D: grid size must be strictly positive'
  endif
  do iid=1,lid-1
    call mpi_send(lhorzn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
    call mpi_send(lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
  end do
else                 !workers
  call mpi_recv(lhorzn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!Everyone must allocate space for the grid of input data
allocate(zn(lzn))    !these are module-scope variables
if (flagcart) then
  allocate(rhon(1))  !not used in Cartesian code so just set to something
  allocate(yn(lhorzn))
  lyn=lhorzn
else
  allocate(rhon(lhorzn))
  allocate(yn(1))    !not used in the axisymmetric code so just initialize to something
  lrhon=lhorzn
end if


!Note that the second dimension ("longitude") is singleton so that we are able to also use these vars for 3D input
allocate(dnO(lzn,1,lhorzn),dnN2(lzn,1,lhorzn),dnO2(lzn,1,lhorzn),dvnrho(lzn,1,lhorzn),dvnz(lzn,1,lhorzn),dTn(lzn,1,lhorzn))


!Define a grid (input data) by assuming that the spacing is constant
if (flagcart) then     !Cartesian neutral simulation
  yn=[ ((real(ihorzn,8)-1._wp)*dhorzn, ihorzn=1,lhorzn) ]
  meanyn=sum(yn,1)/size(yn,1)
  yn=yn-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
else
  rhon=[ ((real(ihorzn,8)-1._wp)*dhorzn, ihorzn=1,lhorzn) ]
end if
zn=[ ((real(izn,8)-1._wp)*cfg%dzn, izn=1,lzn) ]

if (myid==0) then
  if (flagcart) then
    print *, 'Creating neutral grid with y,z extent:',minval(yn),maxval(yn),minval(zn),maxval(zn)
  else
    print *, 'Creating neutral grid with rho,z extent:  ',minval(rhon),maxval(rhon),minval(zn),maxval(zn)
  end if
end if


!Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
phi1=cfg%sourcemlon*pi/180d0
theta1=pi/2d0-cfg%sourcemlat*pi/180d0


!Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
if (myid==0) then
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


      !COMPUTE DISTANCES
      gammarads=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)     !this is actually cos(gamma)
      if (gammarads>1._wp) then     !handles weird precision issues in 2D
        gammarads=1._wp
      else if (gammarads<-1._wp) then
        gammarads=-1._wp
      end if
      gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
      rhoimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''

      !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
      theta3=theta2
      phi3=phi1
      gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
      if (gamma1>1._wp) then     !handles weird precision issues in 2D
        gamma1=1._wp
      else if (gamma1<-1._wp) then
        gamma1=-1._wp
      end if
      gamma1=acos(gamma1)

      gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
      if (gamma2>1._wp) then     !handles weird precision issues in 2D
        gamma2=1._wp
      else if (gamma2<-1._wp) then
        gamma2=-1._wp
      end if
      gamma2=acos(gamma2)

      xp=Re*gamma1
      yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement


      !COMPUTE COORDIANTES FROM DISTANCES
      if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
        yp = -yp            !do we want an abs here to be safe
      end if
      if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than soure point (phi3=phi1)
        xp = -xp
      end if
      phip=atan2(yp,xp)


      if(flagcart) then
        yimat(ix1,ix2,ix3)=yp
      end if


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

      if (flagcart) then
        eyp=-1._wp*x%etheta(ix1,ix2,ix3,:)

        tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e1(ix1,ix2,ix3)=tmpsca

        tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e2(ix1,ix2,ix3)=tmpsca

        tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_eyp_e3(ix1,ix2,ix3)=tmpsca
      else
        erhop=cos(phip)*x%e3(ix1,ix2,ix3,:)+(-1._wp)*sin(phip)*x%etheta(ix1,ix2,ix3,:)     !unit vector for azimuth (referenced from epicenter - not geocenter!!!) in cartesian geocentric-geomagnetic coords.

        tmpvec=erhop*x%e1(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e1(ix1,ix2,ix3)=tmpsca

        tmpvec=erhop*x%e2(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e2(ix1,ix2,ix3)=tmpsca

        tmpvec=erhop*x%e3(ix1,ix2,ix3,:)
        tmpsca=sum(tmpvec)
        proj_erhop_e3(ix1,ix2,ix3)=tmpsca
      end if
    end do
  end do
end do


!Assign values for flat lists of grid points
zi=pack(zimat,.true.)     !create a flat list of grid points to be used by interpolation ffunctions
if (flagcart) then
  yi=pack(yimat,.true.)
else
  rhoi=pack(rhoimat,.true.)
end if


!GRID UNIT VECTORS NO LONGER NEEDED ONCE PROJECTIONS ARE CALCULATED...
call clear_unitvecs(x)


!PRINT OUT SOME BASIC INFO ABOUT THE GRID THAT WE'VE LOADED
if (myid==0 .and. debug) then
  if (flagcart) then
    print *, 'Min/max yn,zn values',minval(yn),maxval(yn),minval(zn),maxval(zn)
    print *, 'Min/max yi,zi values',minval(yi),maxval(yi),minval(zi),maxval(zi)
  else
    print *, 'Min/max rhon,zn values',minval(rhon),maxval(rhon),minval(zn),maxval(zn)
    print *, 'Min/max rhoi,zi values',minval(rhoi),maxval(rhoi),minval(zi),maxval(zi)
  end if

  print *, 'Source lat/long:  ',cfg%sourcemlat,cfg%sourcemlon
  print *, 'Plasma grid lat range:  ',minval(x%glat(:,:,:)),maxval(x%glat(:,:,:))
  print *, 'Plasma grid lon range:  ',minval(x%glon(:,:,:)),maxval(x%glon(:,:,:))
end if

end subroutine gridproj_dneu2D


subroutine gridproj_dneu3D(cfg,x)

!Read in the grid for the neutral data and project unit vectors into the appropriiate directions.
!Also allocate module-scope variables for storing neutral perturbations read in from input files.

type(gemini_cfg), intent(in) :: cfg
type(curvmesh), intent(inout) :: x           !inout to allow deallocation of unit vectors once we are done with them, should consider exporting this to another functino to be called from main program to avoid having x writeable...

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


!Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
phi1=cfg%sourcemlon*pi/180._wp
theta1=pi/2._wp-cfg%sourcemlat*pi/180._wp


!Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
if (myid==0) then
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
      !if (gammarads>1._wp) then     !handles weird precision issues in 2D
      !  gammarads=1._wp
      !else if (gammarads<-1._wp) then
      !  gammarads=-1._wp
      !end if
      !gammarads=acos(gammarads)                     !angle between source location annd field point (in radians)
      !rhoimat(ix1,ix2,ix3)=Re*gammarads    !rho here interpreted as the arc-length defined by angle between epicenter and ``field point''
      !! ZZZ end possibly superfluous block of code...

      !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
      theta3=theta2
      phi3=phi1
      gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
      if (gamma1>1._wp) then     !handles weird precision issues in 2D
        gamma1=1._wp
      else if (gamma1<-1._wp) then
        gamma1=-1._wp
      end if
      gamma1=acos(gamma1)

      gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
      if (gamma2>1._wp) then     !handles weird precision issues in 2D
        gamma2=1._wp
      else if (gamma2<-1._wp) then
        gamma2=-1._wp
      end if
      gamma2=acos(gamma2)
      xp=Re*gamma1
      yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement


      !COMPUTE COORDIANTES FROM DISTANCES
      if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
        yp=-1._wp*yp            !do we want an abs here to be safe
      end if
      if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than soure point (phi3=phi1)
        xp=-1._wp*xp
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

      eyp=-1._wp*x%etheta(ix1,ix2,ix3,:)

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


!GRID UNIT VECTORS NO LONGER NEEDED ONCE PROJECTIONS ARE CALCULATED, so go ahead and free some space
if (myid==0) then
  print*, '...Clearing out unit vectors (after projections)...'
end if
call clear_unitvecs(x)


if(myid==0) then
  print*, 'Projection checking:  ',minval(proj_exp_e1),maxval(proj_exp_e1),minval(proj_exp_e2),maxval(proj_exp_e2), &
                                   minval(proj_exp_e3),maxval(proj_exp_e3)
end if


!Establish the size of the grid based on input file and distribute to workers
if (myid==0) then    !root
  print '(A,/,A)', 'READ neutral size from:', cfg%sourcedir

  call get_simsize3(cfg%sourcedir, lx1=lxnall, lx2all=lynall, lx3all=lzn)

  print *, 'Neutral data has lx,ly,lz size:  ',lxnall,lynall,lzn,' with spacing dx,dy,dz',cfg%dxn,cfg%drhon,cfg%dzn
  if (lxnall < 1 .or. lynall < 1 .or. lzn < 1) then
    write(stderr,*) 'ERROR: reading ' // cfg%sourcedir
    error stop 'neutral:gridproj_dneu3D: grid size must be strictly positive'
  endif

  !root must allocate space for the entire grid of input data - this might be doable one parameter at a time???
  allocate(zn(lzn))        !the z coordinate is never split up in message passing - want to use full altitude range...
  allocate(rhon(1))        !not used in Cartesian or 3D code so just set to something; could likely be left unallocated
  allocate(xnall(lxnall))
  allocate(ynall(lynall))
  allocate(dnOall(lzn,lxnall,lynall),dnN2all(lzn,lxnall,lynall),dnO2all(lzn,lxnall,lynall),dvnrhoall(lzn,lxnall,lynall), &
             dvnzall(lzn,lxnall,lynall),dvnxall(lzn,lxnall,lynall),dTnall(lzn,lxnall,lynall))    !ZZZ - note that these might be deallocated after each read to clean up memory management a bit...


  !calculate the z grid (same for all) and distribute to workers so we can figure out their x-y slabs
  print*, '...creating vertical grid and sending to workers...'
  zn=[ ((real(izn,8)-1._wp)*cfg%dzn, izn=1,lzn) ]    !root calculates and distributes but this is the same for all workers - assmes that the max neutral grid extent in altitude is always less than the plasma grid (should almost always be true)
  maxzn=maxval(zn)
  do iid=1,lid-1
    call mpi_send(lzn,1,MPI_INTEGER,iid,tag%lz,MPI_COMM_WORLD,ierr)
    call mpi_send(zn,lzn,mpi_realprec,iid,tag%zn,MPI_COMM_WORLD,ierr)
  end do


  !Define a global neutral grid (input data) by assuming that the spacing is constant
  ynall=[ ((real(iyn,8)-1._wp)*cfg%drhon, iyn=1,lynall) ]
  meanyn=sum(ynall,1)/size(ynall,1)
  ynall=ynall-meanyn     !the neutral grid should be centered on zero for a cartesian interpolation
  xnall=[ ((real(ixn,8)-1._wp)*cfg%dxn, ixn=1,lxnall) ]
  meanxn=sum(xnall,1)/size(xnall,1)
  xnall=xnall-meanxn     !the neutral grid should be centered on zero for a cartesian interpolation
  print *, 'Created full neutral grid with y,z extent:',minval(xnall),maxval(xnall),minval(ynall), &
                maxval(ynall),minval(zn),maxval(zn)


  !calculate the extent of my piece of the grid using max altitude specified for the neutral grid
  call slabrange(maxzn,ximat,yimat,zimat,cfg%sourcemlat,xnrange,ynrange)
  allocate(extents(0:lid-1,6),indx(0:lid-1,6),slabsizes(0:lid-1,2))
  extents(0,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]


  !receive extents of each of the other workers: extents(lid,6)
  print*, 'Receiving xn and yn ranges from workers...'
  do iid=1,lid-1
    call mpi_recv(xnrange,2,mpi_realprec,iid,tag%xnrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call mpi_recv(ynrange,2,mpi_realprec,iid,tag%ynrange,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    extents(iid,1:6)=[0._wp,maxzn,xnrange(1),xnrange(2),ynrange(1),ynrange(2)]     !need to store values as xnrange overwritten for each worker
    print*, 'Subgrid extents:  ',iid,extents(iid,:)
  end do


  !find index into into neutral arrays for each worker:  indx(lid,6)
  print*, 'Root grid check:  ',ynall(1),ynall(lynall)
  print*, 'Converting ranges to indices...'
  do iid=0,lid-1
    call range2inds(extents(iid,1:6),zn,xnall,ynall,indices)
    indx(iid,1:6)=indices
    print*, 'Subgrid indices',iid,indx(iid,:)
  end do


  !send each worker the sizes for their particular chunk (all different) and send worker that grid chunk
  print*,'Sending sizes and xn,yn subgrids to workers...'
  do iid=1,lid-1
    lxn=indx(iid,4)-indx(iid,3)+1
    lyn=indx(iid,6)-indx(iid,5)+1
    slabsizes(iid,1:2)=[lxn,lyn]
    call mpi_send(lyn,1,MPI_INTEGER,iid,tag%lrho,MPI_COMM_WORLD,ierr)
    call mpi_send(lxn,1,MPI_INTEGER,iid,tag%lx,MPI_COMM_WORLD,ierr)
    allocate(xn(lxn),yn(lyn))
    xn=xnall(indx(iid,3):indx(iid,4))
    yn=ynall(indx(iid,5):indx(iid,6))
    call mpi_send(xn,lxn,mpi_realprec,iid,tag%xn,MPI_COMM_WORLD,ierr)
    call mpi_send(yn,lyn,mpi_realprec,iid,tag%yn,MPI_COMM_WORLD,ierr)
    deallocate(xn,yn)
  end do


  !have root store its part to the full neutral grid
  print*, 'Root is picking out its own subgrid...'
  lxn=indx(0,4)-indx(0,3)+1
  lyn=indx(0,6)-indx(0,5)+1
  slabsizes(0,1:2)=[lxn,lyn]
  allocate(xn(lxn),yn(lyn))
  xn=xnall(indx(0,3):indx(0,4))
  yn=ynall(indx(0,5):indx(0,6))
else                 !workers
  !get teh z-grid from root so we know what the max altitude we have to deal with will be
  call mpi_recv(lzn,1,MPI_INTEGER,0,tag%lz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  allocate(zn(lzn))
  call mpi_recv(zn,lzn,mpi_realprec,0,tag%zn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  maxzn=maxval(zn)


  !calculate the extent of my grid
  call slabrange(maxzn,ximat,yimat,zimat,cfg%sourcemlat,xnrange,ynrange)


  !send ranges to root
  call mpi_send(xnrange,2,mpi_realprec,0,tag%xnrange,MPI_COMM_WORLD,ierr)
  call mpi_send(ynrange,2,mpi_realprec,0,tag%ynrange,MPI_COMM_WORLD,ierr)


  !receive my sizes from root, allocate then receive my pieces of the grid
  call mpi_recv(lxn,1,MPI_INTEGER,0,tag%lx,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(lyn,1,MPI_INTEGER,0,tag%lrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  allocate(xn(lxn),yn(lyn))
  call mpi_recv(xn,lxn,mpi_realprec,0,tag%xn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(yn,lyn,mpi_realprec,0,tag%yn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!AT THIS POINT WE CAN ALLOCATE THE SUBGRID SIZES
allocate(dnO(lzn,lxn,lyn),dnN2(lzn,lxn,lyn),dnO2(lzn,lxn,lyn),dvnrho(lzn,lxn,lyn), &
           dvnz(lzn,lxn,lyn),dvnx(lzn,lxn,lyn),dTn(lzn,lxn,lyn))


!PRINT OUT SOME BASIC INFO ABOUT THE GRID THAT WE'VE LOADED
if (debug) then
  print *, 'Min/max zn,xn,yn values',myid,minval(zn),maxval(zn),minval(xn),maxval(xn),minval(yn),maxval(yn)
  print *, 'Min/max zi,xi,yi values',myid,minval(zi),maxval(zi),minval(xi),maxval(xi),minval(yi),maxval(yi)
  !print *, 'Source lat/long:  ',myid,meanlat,meanlong
  !print *, 'Plasma grid lat range:  ',myid,minval(x%glat(:,:,:)),maxval(x%glat(:,:,:))
  !print *, 'Plasma grid lon range:  ',myid,minval(x%glon(:,:,:)),maxval(x%glon(:,:,:))
end if

end subroutine gridproj_dneu3D


subroutine read_dneu2D(tprev,tnext,t,dtneu,dt,neudir,ymdtmp,UTsectmp,flagcart)

! This subroutine reads in neutral frame data, if necessary and puts the results
! in a module-scope variable for later use
! This version is for 2D source neutral data

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt        !times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir                    !directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp          !storage space for incrementing date without overwriting ymdnext...
real(wp), intent(out) :: UTsectmp
logical, intent(in) :: flagcart

integer :: iid,ierr
integer :: lhorzn                        !number of horizontal grid points


if (flagcart) then
  lhorzn=lyn
else
  lhorzn=lrhon
end if


if (myid==0) then    !root
  !read in the data from file
  if(debug) print *, 'neutral.f90:read_dneu2D: tprev,tnow,tnext:  ',tprev,t+dt/2,tnext
  ymdtmp=ymdnext
  UTsectmp=UTsecnext
  call dateinc(dtneu,ymdtmp,UTsectmp)    !get the date for "next" params

  call get_neutral2(date_filename(neudir,ymdtmp,UTsectmp), &
    dnO,dnN2,dnO2,dvnrho,dvnz,dTn)

  if (debug) then
    print *, 'Min/max values for dnO:  ',minval(dnO),maxval(dnO)
    print *, 'Min/max values for dnN:  ',minval(dnN2),maxval(dnN2)
    print *, 'Min/max values for dnO:  ',minval(dnO2),maxval(dnO2)
    print *, 'Min/max values for dvnrho:  ',minval(dvnrho),maxval(dvnrho)
    print *, 'Min/max values for dvnz:  ',minval(dvnz),maxval(dvnz)
    print *, 'Min/max values for dTn:  ',minval(dTn),maxval(dTn)
  endif

  if (.not. all(ieee_is_finite(dnO))) error stop 'dnO: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnN2))) error stop 'dnN2: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnO2))) error stop 'dnO2: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnrho))) error stop 'dvnrho: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnz))) error stop 'dvnz: non-finite value(s)'
  if (.not. all(ieee_is_finite(dTn))) error stop 'dTn: non-finite value(s)'

  !send a full copy of the data to all of the workers
  do iid=1,lid-1
    call mpi_send(dnO,lhorzn*lzn,mpi_realprec,iid,tag%dnO,MPI_COMM_WORLD,ierr)
    call mpi_send(dnN2,lhorzn*lzn,mpi_realprec,iid,tag%dnN2,MPI_COMM_WORLD,ierr)
    call mpi_send(dnO2,lhorzn*lzn,mpi_realprec,iid,tag%dnO2,MPI_COMM_WORLD,ierr)
    call mpi_send(dTn,lhorzn*lzn,mpi_realprec,iid,tag%dTn,MPI_COMM_WORLD,ierr)
    call mpi_send(dvnrho,lhorzn*lzn,mpi_realprec,iid,tag%dvnrho,MPI_COMM_WORLD,ierr)
    call mpi_send(dvnz,lhorzn*lzn,mpi_realprec,iid,tag%dvnz,MPI_COMM_WORLD,ierr)
  end do
else     !workers
  !receive a full copy of the data from root
  call mpi_recv(dnO,lhorzn*lzn,mpi_realprec,0,tag%dnO,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnN2,lhorzn*lzn,mpi_realprec,0,tag%dnN2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnO2,lhorzn*lzn,mpi_realprec,0,tag%dnO2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dTn,lhorzn*lzn,mpi_realprec,0,tag%dTn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnrho,lhorzn*lzn,mpi_realprec,0,tag%dvnrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnz,lhorzn*lzn,mpi_realprec,0,tag%dvnz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!DO SPATIAL INTERPOLATION OF EACH PARAMETER (COULD CONSERVE SOME MEMORY BY NOT STORING DVNRHOIPREV AND DVNRHOINEXT, ETC.)
if (myid==lid/2 .and. debug) then
  print*, 'neutral data size:  ',lhorzn,lzn,lid
  print *, 'Min/max values for dnO:  ',minval(dnO),maxval(dnO)
  print *, 'Min/max values for dnN:  ',minval(dnN2),maxval(dnN2)
  print *, 'Min/max values for dnO:  ',minval(dnO2),maxval(dnO2)
  print *, 'Min/max values for dvnrho:  ',minval(dvnrho),maxval(dvnrho)
  print *, 'Min/max values for dvnz:  ',minval(dvnz),maxval(dvnz)
  print *, 'Min/max values for dTn:  ',minval(dTn),maxval(dTn)
  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
end if

end subroutine read_dneu2D


subroutine read_dneu3D(tprev,tnext,t,dtneu,dt,neudir,ymdtmp,UTsectmp)

! This subroutine reads in neutral frame data, if necessary and puts the results
! in a module-scope variable for later use
! this version is for 3D neutral source data

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt        !times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir                    !directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp          !storage space for incrementing date without overwriting ymdnext...
real(wp), intent(out) :: UTsectmp

integer :: iid,ierr
integer :: lhorzn                        !number of horizontal grid points
real(wp), dimension(:,:,:), allocatable :: parmtmp    !temporary resizable space for subgrid neutral data


lhorzn=lyn


if (myid==0) then    !root
  !read in the data from file
  if(debug) print *, 'tprev,tnow,tnext:  ',tprev,t+dt/2d0,tnext
  ymdtmp=ymdnext
  UTsectmp=UTsecnext
  call dateinc(dtneu,ymdtmp,UTsectmp)                !get the date for "next" params

  call get_neutral3(date_filename(neudir,ymdtmp,UTsectmp), &
    dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)

  if (debug) then
    print *, 'Min/max values for dnOall:  ',minval(dnOall),maxval(dnOall)
    print *, 'Min/max values for dnNall:  ',minval(dnN2all),maxval(dnN2all)
    print *, 'Min/max values for dnO2all:  ',minval(dnO2all),maxval(dnO2all)
    print *, 'Min/max values for dvnxall:  ',minval(dvnxall),maxval(dvnxall)
    print *, 'Min/max values for dvnrhoall:  ',minval(dvnrhoall),maxval(dvnrhoall)
    print *, 'Min/max values for dvnzall:  ',minval(dvnzall),maxval(dvnzall)
    print *, 'Min/max values for dTnall:  ',minval(dTnall),maxval(dTnall)
  endif

  if (.not. all(ieee_is_finite(dnOall))) error stop 'dnOall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnN2all))) error stop 'dnN2all: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnO2all))) error stop 'dnO2all: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnxall))) error stop 'dvnxall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnrhoall))) error stop 'dvnrhoall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnzall))) error stop 'dvnzall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dTnall))) error stop 'dTnall: non-finite value(s)'


  !in the 3D case we cannnot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
  do iid=1,lid-1
    allocate(parmtmp(lzn,slabsizes(iid,1),slabsizes(iid,2)))    !get space for the parameter for this worker

    parmtmp=dnOall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnO,MPI_COMM_WORLD,ierr)

    parmtmp=dnN2all(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnN2,MPI_COMM_WORLD,ierr)

    parmtmp=dnO2all(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnO2,MPI_COMM_WORLD,ierr)

    parmtmp=dTnall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dTn,MPI_COMM_WORLD,ierr)

    parmtmp=dvnrhoall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnrho,MPI_COMM_WORLD,ierr)

    parmtmp=dvnzall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnz,MPI_COMM_WORLD,ierr)

    parmtmp=dvnxall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnx,MPI_COMM_WORLD,ierr)

    deallocate(parmtmp)
  end do

  !root needs to grab its piece of data
  dnO=dnOall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dnN2=dnN2all(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dnO2=dnO2all(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dTn=dTnall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnrho=dvnrhoall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnz=dvnzall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnx=dvnxall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
else     !workers
  !receive a subgrid copy of the data from root
  call mpi_recv(dnO,lzn*lxn*lyn,mpi_realprec,0,tag%dnO,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnN2,lzn*lxn*lyn,mpi_realprec,0,tag%dnN2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnO2,lzn*lxn*lyn,mpi_realprec,0,tag%dnO2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dTn,lzn*lxn*lyn,mpi_realprec,0,tag%dTn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnrho,lzn*lxn*lyn,mpi_realprec,0,tag%dvnrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnz,lzn*lxn*lyn,mpi_realprec,0,tag%dvnz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnx,lzn*lxn*lyn,mpi_realprec,0,tag%dvnx,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


if (myid==lid/2 .and. debug) then
  print*, 'neutral data size:  ',myid,lzn,lxn,lyn
  print *, 'Min/max values for dnO:  ',myid,minval(dnO),maxval(dnO)
  print *, 'Min/max values for dnN:  ',myid,minval(dnN2),maxval(dnN2)
  print *, 'Min/max values for dnO2:  ',myid,minval(dnO2),maxval(dnO2)
  print *, 'Min/max values for dvnx:  ',myid,minval(dvnx),maxval(dvnx)
  print *, 'Min/max values for dvnrho:  ',myid,minval(dvnrho),maxval(dvnrho)
  print *, 'Min/max values for dvnz:  ',myid,minval(dvnz),maxval(dvnz)
  print *, 'Min/max values for dTn:  ',myid,minval(dTn),maxval(dTn)
!  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
end if

end subroutine read_dneu3D


subroutine spaceinterp_dneu2D(flagcart)

!must take into account the type of interpolation that is being done
logical, intent(in) :: flagcart

integer :: lhorzn
real(wp), dimension(:,:), allocatable :: tmpinterp
real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data


!Array for packing neutral data
if (flagcart) then
  lhorzn=lyn
else
  lhorzn=lrhon
end if
allocate(tmpinterp(lzn,lhorzn))


if(flagcart) then
  tmpinterp=dnO(:,1,:)                    !pack into 2D array for interp2
  parami=interp2(zn,yn,tmpinterp,zi,yi)         !interp to temp var.
  dnOiprev=dnOinext                       !save new previous
  dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

  tmpinterp=dnN2(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dnN2iprev=dnN2inext
  dnN2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dnO2(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dnO2iprev=dnO2inext
  dnO2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnrho(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
  dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnz(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dvnziprev=dvnzinext
  dvnzinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dTn(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dTniprev=dTninext
  dTninext=reshape(parami,[lx1,lx2,lx3])
else
  tmpinterp=dnO(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)     !interp to temp var.
  dnOiprev=dnOinext                       !save new previous
  dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

  tmpinterp=dnN2(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dnN2iprev=dnN2inext
  dnN2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dnO2(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dnO2iprev=dnO2inext
  dnO2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnrho(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dvnrhoiprev=dvnrhoinext
  dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnz(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dvnziprev=dvnzinext
  dvnzinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dTn(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dTniprev=dTninext
  dTninext=reshape(parami,[lx1,lx2,lx3])
end if

!MORE DIAG
if (myid==lid/2) then
  print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvrhoi:  ',minval(dvnrhoinext),maxval(dvnrhoinext)
  print *, 'Min/max values for dvnzi:  ',minval(dvnzinext),maxval(dvnzinext)
  print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
end if


!ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
!structure)
dvn1iprev=dvn1inext   !save the old data
dvn2iprev=dvn2inext
dvn3iprev=dvn3inext
if(flagcart) then
  dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
  dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2
  dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3
else
  dvn1inext=dvnrhoinext*proj_erhop_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates
  dvn2inext=dvnrhoinext*proj_erhop_e2+dvnzinext*proj_ezp_e2
  dvn3inext=dvnrhoinext*proj_erhop_e3+dvnzinext*proj_ezp_e3
end if

!MORE DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvn1i:  ',minval(dvn1inext),maxval(dvn1inext)
  print *, 'Min/max values for dvn2i:  ',minval(dvn2inext),maxval(dvn2inext)
  print *, 'Min/max values for dvn3i:  ',minval(dvn3inext),maxval(dvn3inext)
  print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
end if


!CLEAR ALLOCATED VARS
deallocate(tmpinterp)

end subroutine spaceinterp_dneu2D


subroutine spaceinterp_dneu3D()

!performs spatial interpolation for 3D input neutral data from MAGIC or some other source

real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data


!INTERPOLATE IN THREE DIMENSIONS
parami=interp3(zn,xn,yn,dnO,zi,xi,yi)         !interp to temp var.
dnOiprev=dnOinext                       !save new previous
dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

parami=interp3(zn,xn,yn,dnN2,zi,xi,yi)
dnN2iprev=dnN2inext
dnN2inext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dnO2,zi,xi,yi)
dnO2iprev=dnO2inext
dnO2inext=reshape(parami,[lx1,lx2,lx3])

!ZZZ - do we want to make dvnrho-->dvny???
parami=interp3(zn,xn,yn,dvnrho,zi,xi,yi)
dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnz,zi,xi,yi)
dvnziprev=dvnzinext
dvnzinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnx,zi,xi,yi)
dvnxiprev=dvnxinext
dvnxinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dTn,zi,xi,yi)
dTniprev=dTninext
dTninext=reshape(parami,[lx1,lx2,lx3])


!MORE DIAG
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dnOi:  ',myid,minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',myid,minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',myid,minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvrhoi:  ',myid,minval(dvnrhoinext),maxval(dvnrhoinext)
  print *, 'Min/max values for dvnzi:  ',myid,minval(dvnzinext),maxval(dvnzinext)
  print *, 'Min/max values for dvnxi:  ',myid,minval(dvnxinext),maxval(dvnxinext)
  print *, 'Min/max values for dTni:  ',myid,minval(dTninext),maxval(dTninext)
end if


!ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
!structure)
dvn1iprev=dvn1inext   !save the old data for the rotated vectors
dvn2iprev=dvn2inext
dvn3iprev=dvn3inext
dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1+dvnxinext*proj_exp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2+dvnxinext*proj_exp_e2
dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3+dvnxinext*proj_exp_e3


!MORE DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dvn1i:  ',myid,minval(dvn1inext),maxval(dvn1inext)
  print *, 'Min/max values for dvn2i:  ',myid,minval(dvn2inext),maxval(dvn2inext)
  print *, 'Min/max values for dvn3i:  ',myid,minval(dvn3inext),maxval(dvn3inext)
end if

end subroutine spaceinterp_dneu3D


subroutine timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!interpolation in time - no sensitive to dimensionality of the input neutral data so this can be
! the same for 2D vs. 3D

real(wp), intent(in) :: t,dt
real(wp), dimension(:,:,:), intent(out) :: dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow

integer :: ix1,ix2,ix3
real(wp) :: slope

do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      slope=(dnOinext(ix1,ix2,ix3)-dnOiprev(ix1,ix2,ix3))/(tnext-tprev)
      dnOinow(ix1,ix2,ix3)=dnOiprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnN2inext(ix1,ix2,ix3)-dnN2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnN2inow(ix1,ix2,ix3)=dnN2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnO2inext(ix1,ix2,ix3)-dnO2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnO2inow(ix1,ix2,ix3)=dnO2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn1inext(ix1,ix2,ix3)-dvn1iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn1inow(ix1,ix2,ix3)=dvn1iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn2inext(ix1,ix2,ix3)-dvn2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn2inow(ix1,ix2,ix3)=dvn2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn3inext(ix1,ix2,ix3)-dvn3iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn3inow(ix1,ix2,ix3)=dvn3iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dTninext(ix1,ix2,ix3)-dTniprev(ix1,ix2,ix3))/(tnext-tprev)
      dTninow(ix1,ix2,ix3)=dTniprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)
    end do
  end do
end do


!SOME BASIC DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'tprev,t,tnext:  ',myid,tprev,t+dt/2d0,tnext
  print *, 'Min/max values for dnOinow:  ',myid,minval(dnOinow),maxval(dnOinow)
  print *, 'Min/max values for dnN2inow:  ',myid,minval(dnN2inow),maxval(dnN2inow)
  print *, 'Min/max values for dnO2inow:  ',myid,minval(dnO2inow),maxval(dnO2inow)
  print *, 'Min/max values for dvn1inow:  ',myid,minval(dvn1inow),maxval(dvn1inow)
  print *, 'Min/max values for dvn2inow:  ',myid,minval(dvn2inow),maxval(dvn2inow)
  print *, 'Min/max values for dvn3inow:  ',myid,minval(dvn3inow),maxval(dvn3inow)
  print *, 'Min/max values for dTninow:  ',myid,minval(dTninow),maxval(dTninow)
end if

end subroutine timeinterp_dneu


subroutine neutral_update(nn,Tn,vn1,vn2,vn3)

! adds stored base and perturbation neutral atmospheric parameters (these are module-scope parameters so not needed as input)

real(wp), dimension(:,:,:,:), intent(out) :: nn
real(wp), dimension(:,:,:), intent(out) :: Tn
real(wp), dimension(:,:,:), intent(out) :: vn1,vn2,vn3


!> background neutral parameters
nn=nnmsis
Tn=Tnmsis
vn1=vn1base
vn2=vn2base
vn3=vn3base

!> perturbations, if used
if (allocated(zn)) then
  nn(:,:,:,1)=nn(:,:,:,1)+dnOinow
  nn(:,:,:,2)=nn(:,:,:,2)+dnN2inow
  nn(:,:,:,3)=nn(:,:,:,3)+dnO2inow
  nn(:,:,:,1)=max(nn(:,:,:,1),1._wp)
  nn(:,:,:,2)=max(nn(:,:,:,2),1._wp)
  nn(:,:,:,3)=max(nn(:,:,:,3),1._wp)
  !! note we are not adjust derived densities like NO since it's not clear how they may be related to major
  !! species perturbations.

  Tn=Tn+dTninow
  Tn=max(Tn,50._wp)

  vn1=vn1+dvn1inow
  vn2=vn2+dvn2inow
  vn3=vn3+dvn3inow
end if

end subroutine neutral_update


subroutine slabrange(maxzn,ximat,yimat,zimat,sourcemlat,xnrange,ynrange)

!takes in a subgrid and the max altitude of interest for neutral interpolation and then computes
!what the maximum xn and yn will be for that slab
! ZZZ - also this is specific to dipole grids right now...

real(wp), intent(in) :: maxzn
real(wp), dimension(:,:,:), intent(in) :: ximat,yimat,zimat
real(wp), intent(in) :: sourcemlat
real(wp), dimension(2), intent(out) :: xnrange,ynrange     !for min and max

real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
integer :: lx1tmp
integer, dimension(size(ximat,2),size(ximat,3)) :: ix1stmp
integer :: ix1tmp
logical :: flagSH
integer :: ix1


!in what hemisphere is our source?
if (sourcemlat<=0) then
  flagSH=.true.
else
  flagSH=.false.
end if


!peel the grid in half (source hemisphere if closed dipole)
if (gridflag==0) then    !closed dipole grid

  ix1=maxloc(pack(zimat(:,1,1),.true.),1)    !apex is by definition the highest altitude along a given field line
  if (flagSH) then
    lx1tmp=ix1                  !first piece of arrays
  else
    lx1tmp=lx1-ix1    !second (end) piece of arrays
  end if
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers???

  if(flagSH) then    !southern hemisphere
    xitmp=ximat(1:ix1,1:lx2,1:lx3)          !select beginning of the array - the southern half
    yitmp=yimat(1:ix1,1:lx2,1:lx3)
    zitmp=zimat(1:ix1,1:lx2,1:lx3)
  else               !northern hemisphere
      xitmp=ximat(ix1+1:lx1,1:lx2,1:lx3)    !select end half of the array
      yitmp=yimat(ix1+1:lx1,1:lx2,1:lx3)
      zitmp=zimat(ix1+1:lx1,1:lx2,1:lx3)
  end if
else     !this is not an interhemispheric grid so our approach is to just use all of the data
  lx1tmp=lx1
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers?
  xitmp=ximat(1:lx1,1:lx2,1:lx3)
  yitmp=yimat(1:lx1,1:lx2,1:lx3)
  zitmp=zimat(1:lx1,1:lx2,1:lx3)
!  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
end if


!the min and max x are simply determined by longitude...
xnrange(1)=minval(xitmp)
xnrange(2)=maxval(xitmp)


!situation is more complicated for latitude due to dipole grid, need to determine by L-shell
if (flagSH) then
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0, just use the first longitude slice since they will all have the same L-shell-field line relations
  ynrange(2)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(1)=yitmp(ix1,lx2,1)
else    !things are swapped around in NH
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0; this is the southernmost edge of the neutral slab we need
  ynrange(1)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(2)=yitmp(ix1,lx2,1)
end if

deallocate(xitmp,yitmp,zitmp)

end subroutine slabrange


subroutine range2inds(ranges,zn,xnall,ynall,indices)

!determine where the slab described by ranges falls within the global neutral grid

real(wp), dimension(6), intent(in) :: ranges
real(wp), dimension(:), intent(in) :: zn,xnall,ynall
integer, dimension(6), intent(out) :: indices

real(wp) :: minzn,maxzn,minxn,maxxn,minyn,maxyn
integer :: ixn,iyn


!for clarity
minzn=ranges(1)
maxzn=ranges(2)
minxn=ranges(3)
maxxn=ranges(4)
minyn=ranges(5)
maxyn=ranges(6)


!always use the full z-range
indices(1)=1
indices(2)=lzn


!x-range
ixn=1
do while (ixn<lxnall .and. xnall(ixn)<minxn)
  ixn=ixn+1
end do
indices(3)=max(ixn-1,1)    !just to be sure go back one index so that we cover the min range, don't let index fall below zero
do while (ixn<lxnall .and. xnall(ixn)<maxxn)
  ixn=ixn+1
end do
indices(4)=ixn


!y-range
iyn=1
do while (iyn<lynall .and. ynall(iyn)<minyn)
  iyn=iyn+1
end do
indices(5)=max(iyn-1,1)    !just to be sure go back one index so that we cover the min range
do while (iyn<lynall .and. ynall(iyn)<maxyn)
  iyn=iyn+1
end do
indices(6)=iyn

print*, '!!!!!!!!!!!!!!!!!'
print*, myid
print*, ranges
print*, indices
print*, lxnall,lynall
print*, xnall(indices(3)),xnall(indices(4))
print*, ynall(indices(5)),ynall(indices(6))
print*, '!!!!!!!!!!!!!!!!!'


!corner cases - range is not at all within the neutral grid...  Manifests as both indices being either 1 or lxi, interpolation should zero these out...

end subroutine range2inds


subroutine make_dneu()

!ZZZ - could make this take in type of neutral interpolation and do allocations accordingly

!allocate and compute plasma grid z,rho locations and space to save neutral perturbation variables and projection factors
allocate(zi(lx1*lx2*lx3),rhoi(lx1*lx2*lx3))
allocate(yi(lx1*lx2*lx3))
allocate(xi(lx1*lx2*lx3))
allocate(proj_erhop_e1(lx1,lx2,lx3),proj_ezp_e1(lx1,lx2,lx3),proj_erhop_e2(lx1,lx2,lx3),proj_ezp_e2(lx1,lx2,lx3), &
         proj_erhop_e3(lx1,lx2,lx3),proj_ezp_e3(lx1,lx2,lx3))
allocate(proj_eyp_e1(lx1,lx2,lx3),proj_eyp_e2(lx1,lx2,lx3),proj_eyp_e3(lx1,lx2,lx3))
allocate(proj_exp_e1(lx1,lx2,lx3),proj_exp_e2(lx1,lx2,lx3),proj_exp_e3(lx1,lx2,lx3))
allocate(dnOiprev(lx1,lx2,lx3),dnN2iprev(lx1,lx2,lx3),dnO2iprev(lx1,lx2,lx3),dvnrhoiprev(lx1,lx2,lx3), &
         dvnziprev(lx1,lx2,lx3),dTniprev(lx1,lx2,lx3),dvn1iprev(lx1,lx2,lx3),dvn2iprev(lx1,lx2,lx3), &
         dvn3iprev(lx1,lx2,lx3))
allocate(dvnxiprev(lx1,lx2,lx3))
allocate(dnOinext(lx1,lx2,lx3),dnN2inext(lx1,lx2,lx3),dnO2inext(lx1,lx2,lx3),dvnrhoinext(lx1,lx2,lx3), &
         dvnzinext(lx1,lx2,lx3),dTninext(lx1,lx2,lx3),dvn1inext(lx1,lx2,lx3),dvn2inext(lx1,lx2,lx3), &
         dvn3inext(lx1,lx2,lx3))
allocate(dvnxinext(lx1,lx2,lx3))
allocate(nnmsis(lx1,lx2,lx3,lnchem),Tnmsis(lx1,lx2,lx3),vn1base(lx1,lx2,lx3),vn2base(lx1,lx2,lx3),vn3base(lx1,lx2,lx3))
allocate(dnOinow(lx1,lx2,lx3),dnN2inow(lx1,lx2,lx3),dnO2inow(lx1,lx2,lx3),dvn1inow(lx1,lx2,lx3),dvn2inow(lx1,lx2,lx3), &
           dvn3inow(lx1,lx2,lx3), dTninow(lx1,lx2,lx3))

!start everyone out at zero
zi=0d0; rhoi=0d0; yi=0d0;
xi=0d0
proj_erhop_e1=0d0; proj_ezp_e1=0d0;
proj_erhop_e2=0d0; proj_ezp_e2=0d0;
proj_erhop_e3=0d0; proj_ezp_e3=0d0;
proj_eyp_e1=0d0; proj_eyp_e2=0d0; proj_eyp_e3=0d0;
proj_exp_e1=0d0; proj_exp_e2=0d0; proj_exp_e3=0d0;
dnOiprev=0d0; dnN2iprev=0d0; dnO2iprev=0d0; dTniprev=0d0; dvnrhoiprev=0d0; dvnziprev=0d0;
dvn1iprev=0d0; dvn2iprev=0d0; dvn3iprev=0d0;
dvnxiprev=0d0
dnOinext=0d0; dnN2inext=0d0; dnO2inext=0d0; dTninext=0d0; dvnrhoinext=0d0; dvnzinext=0d0;
dvn1inext=0d0; dvn2inext=0d0; dvn3inext=0d0;
dvnxinext=0d0
nnmsis=0d0; Tnmsis=0d0; vn1base=0d0; vn2base=0d0; vn3base=0d0
dnOinow=0d0; dnN2inow=0d0; dnO2inow=0d0; dTninow=0d0; dvn1inow=0d0; dvn2inow=0d0; dvn3inow=0d0

!now initialize some module variables
tprev=0d0; tnext=0d0

end subroutine make_dneu


subroutine clear_dneu

!stuff allocated at beginning of program
deallocate(zi,rhoi)
deallocate(yi)
deallocate(proj_erhop_e1,proj_ezp_e1,proj_erhop_e2,proj_ezp_e2, &
         proj_erhop_e3,proj_ezp_e3)
deallocate(proj_eyp_e1,proj_eyp_e2,proj_eyp_e3)
deallocate(dnOiprev,dnN2iprev,dnO2iprev,dvnrhoiprev,dvnziprev,dTniprev,dvn1iprev,dvn2iprev,dvn3iprev)
deallocate(dnOinext,dnN2inext,dnO2inext,dvnrhoinext,dvnzinext,dTninext,dvn1inext,dvn2inext,dvn3inext)
deallocate(nnmsis,Tnmsis,vn1base,vn2base,vn3base)
deallocate(dnOinow,dnN2inow,dnO2inow,dTninow,dvn1inow,dvn2inow,dvn3inow)

!check whether any other module variables were allocated and deallocate accordingly
if (allocated(zn) ) then    !if one is allocated, then they all are
  deallocate(zn)
  deallocate(dnO,dnN2,dnO2,dvnrho,dvnz,dTn)
end if
if (allocated(rhon)) then
  deallocate(rhon)
end if
if (allocated(yn)) then
  deallocate(yn)
end if
if (allocated(extents)) then
  deallocate(extents,indx,slabsizes)
end if
if (allocated(dvnx)) then
  deallocate(dvnx)
end if
if (allocated(xn)) then
  deallocate(xn)
end if
if (allocated(xnall)) then
  deallocate(xnall,ynall)
  deallocate(dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)
end if

end subroutine clear_dneu

end module neutral
