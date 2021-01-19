Program MagCalc
!! Need program statement for FORD
!! THIS IS THE MAIN PROGRAM FOR COMPUTING MAGNETIC FIELDS
!! FROM OUTPUT FROM A SIMULATIONS DONE BY GEMINI3D.
!! THIS PROGRAM VERY MUCH MIRRORS THE SETUP OF THE MAIN GEMINI.F90 CODE.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use magcalc_cli, only : cli
use phys_consts, only : pi,mu0, wp, re, debug
use grid, only : lx1, lx2, lx3, read_grid, clear_grid, lx2all,lx3all,grid_size
use mesh, only : curvmesh
use timeutils, only : dateinc,find_time_elapsed
use config, only : gemini_cfg
use io, only : input_plasma_currents,create_outdir_mag,output_magfields
use mpimod, only: mpi_sum, mpi_comm_world, &
mpibreakdown, process_grid_auto, mpi_manualgrid, halo_end, &
mpi_cfg, mpi_realprec, tag=>gemini_mpi
use h5fortran, only : hdf5_file
use pathlib, only : get_suffix

implicit none (type, external)

external :: mpi_reduce

!> VARIABLES READ IN FROM CONFIG FILE

real(wp) :: UTsec
!! UT (s)
integer, dimension(3) :: ymd
!! year, month, day (current)

type(gemini_cfg) :: cfg
!! most user parameters

!! GRID STRUCTURE
type(curvmesh) :: x
!! structure containg grid locations, finite differences, etc.:  see grid module for details

!STATE VARIABLES
real(wp), dimension(:,:,:), allocatable :: J1,J2,J3      !electrodynamic state variables

!TEMPORAL VARIABLES
real(wp) :: t=0, dt      !time from beginning of simulation (s) and time step (s)
real(wp) :: tout    !time for next output and time between outputs
real(wp) :: tstart,tfin   !temp. vars. for measuring performance of code blocks
integer :: it,isp        !time and species loop indices

!WORK ARRAYS
integer :: flag2D
real(wp), dimension(:,:,:), allocatable :: xp,yp,zp     !source coordinates computed from simulation sub-grid
real(wp), dimension(:), allocatable  :: xf,yf,zf        !field point coordinates (flat list)
real(wp), dimension(:), allocatable :: r,theta,phi
real(wp), dimension(:,:,:), allocatable :: dV
integer :: ipoints,lpoints    !size of field point arrays
real(wp), dimension(:,:,:), allocatable :: proj_e1er,proj_e2er,proj_e3er
real(wp), dimension(:,:,:), allocatable :: proj_e1etheta,proj_e2etheta,proj_e3etheta
real(wp), dimension(:,:,:), allocatable :: proj_e1ephi,proj_e2ephi,proj_e3ephi
real(wp), dimension(:,:,:), allocatable :: Jx,Jy,Jz
real(wp), dimension(:,:,:), allocatable :: Rx,Ry,Rz,Rcubed
real(wp), dimension(:,:,:), allocatable :: integrand,integrandavg
real(wp), dimension(:,:,:), allocatable :: alt
real(wp), dimension(:), allocatable :: Br,Btheta,Bphi
real(wp), dimension(:), allocatable :: Brall,Bthetaall,Bphiall
real(wp), dimension(:,:), allocatable :: Jxend,Jyend,Jzend,Rxend,Ryend,Rzend,Rcubedend,dVend
real(wp), dimension(:,:), allocatable :: integrandend
real(wp), dimension(:,:), allocatable :: integrandavgend

real(wp), dimension(:,:), allocatable :: Jxtop,Jytop,Jztop,Rxtop,Rytop,Rztop,Rcubedtop,dVtop
real(wp), dimension(:,:), allocatable :: integrandtop
real(wp), dimension(:,:), allocatable :: integrandavgtop

integer :: ix1,ix2,ix3
real(wp) :: rmean,thetamean

!! FOR SPECIFYING THE PROCESS GRID
integer :: lid2in,lid3in

!! REGULATOR FOR 1/R^3
!real(wp), parameter :: R3min=1d11     !works well for 192x192 in situ
real(wp), parameter :: R3min=1d9
real(wp), parameter :: Rmin=5d3

integer :: ierr

!! for keeping track of start and end times requested by the user
integer, dimension(3) :: ymdstart,ymdend,ymdfinal
real(wp) :: UTsecstart,UTsecend,telend,UTsecfinal
real(wp) :: h1avg,h2avg,h3avg
real(wp), dimension(:,:,:), allocatable :: Rmag


!! --- MAIN PROGRAM

!> get command line parameters and simulation config
call cli(cfg,lid2in,lid3in,debug,ymdstart,UTsecstart,ymdend,UTsecend)


!> set the duration and start time for the magnetic field calculations based on what is give to cli
telend=0;
if (any(ymdstart>0)) then    ! user-specified custom start and end times
  !! count time elapsed up to the requested start time, using output cadence
  UTsec=UTsecstart
  ymd=ymdstart
else   ! start at the beginning of the simulation
  UTsec=cfg%UTsec0
  ymd = cfg%ymd0
end if
if (any(ymdend>0)) then    ! user specified end time (exclusive)
  !! count time elapsed between requested end time and actual simulation end time
  telend=find_time_elapsed(ymd,UTsec,ymdend,UTsecend,cfg%dtout)
else    ! assume user wants to run until the end
  ymdfinal=cfg%ymd0
  UTsecfinal=cfg%UTsec0
  call dateinc(cfg%tdur,ymdfinal,UTsecfinal)
  telend=find_time_elapsed(ymd,UTsec,ymdfinal,UTsecfinal,cfg%dtout)
end if
cfg%tdur=telend


!ESTABLISH A PROCESS GRID
!call grid_size(cfg%indatsize)
!call process_grid_auto(lx2all,lx3all)    !following grid_size these are in scope
!!CHECK THE GRID SIZE AND ESTABLISH A PROCESS GRID
call grid_size(cfg%indatsize)
if (lid2in == -1) then
  !! try to decide the process grid ourself
  call process_grid_auto(lx2all,lx3all)
else
  !! user specified process grid
  call mpi_manualgrid(lx2all,lx3all,lid2in,lid3in)
end if
print '(A, I0, A1, I0)', 'process grid (Number MPI processes) x2, x3:  ',mpi_cfg%lid2, ' ', mpi_cfg%lid3
print '(A, I0, A, I0, A1, I0)', 'Process:',mpi_cfg%myid,' at process grid location: ',mpi_cfg%myid2,' ',mpi_cfg%myid3

!> LOAD UP THE GRID STRUCTURE/MODULE VARS. FOR THIS SIMULATION - THIS ALSO PERMUTES DIMENSIONS OF 2D GRID, IF NEEDED
if (mpi_cfg%myid==0) then
  print*, 'Process grid established; reading in full grid file now...'
end if
call read_grid(cfg%indatsize,cfg%indatgrid, cfg%flagperiodic,x)
!! read in a previously generated grid from filename listed in input file, distribute subgrids to individual workers
if (lx2==1) then
  flag2D=1
  !! a 2D simulation was done, which changes how the integrations go...
else
  flag2D=0
end if

! FIXME:  need to copy the input grid file into the output directory

! Alert the user as to the total simluation time to be computed by magcalc
if (mpi_cfg%myid==0) print*, 'Magcalc total simulation time coverage:  ',telend


!SET UP DIRECTORY TO STORE OUTPUT FILES
if (mpi_cfg%myid==0) then
  call create_outdir_mag(cfg%outdir, cfg%fieldpointfile)
end if


!ALLOCATE ARRAYS (AT THIS POINT ALL SIZES ARE SET FOR EACH PROCESS SUBGRID)
allocate(J1(lx1,lx2,lx3),J2(lx1,lx2,lx3),J3(lx1,lx2,lx3))
allocate(Jx(lx1,lx2,lx3),Jy(lx1,lx2,lx3),Jz(lx1,lx2,lx3))


!NOW DEAL WITH THE UNPRIMED COORDINATES
block
  integer :: u
  type(hdf5_file) :: hf
  character(:), allocatable :: suffix

  suffix = get_suffix(cfg%indatsize)

  select case (suffix)
  case ('.dat')
    open(newunit=u,file=cfg%fieldpointfile,status='old',form='unformatted',access='stream',action='read')
    read(u) lpoints    !size of coordinates for field points
    if (mpi_cfg%myid==0) print *, 'magcalc.f90 --> Number of field points:  ',lpoints
    allocate(r(lpoints),theta(lpoints),phi(lpoints))
    read(u) r,theta,phi
    close(u)
  case ('.h5')
    !! hdf5 file input
    call hf%initialize(cfg%fieldpointfile, status='old', action='r')

    call hf%read('/lpoints',lpoints)
    allocate(r(lpoints), theta(lpoints), phi(lpoints))
    call hf%read('/r',r)
    call hf%read('/theta',theta)
    call hf%read('/phi',phi)
  case default
    error stop 'unrecognized input field point file type'
  end select
end block

if (mpi_cfg%myid==0) print *, 'magcalc.f90 --> Range of r,theta,phi',minval(r),maxval(r),minval(theta), &
                           maxval(theta),minval(phi),maxval(phi)
rmean=sum(r)/size(r)
thetamean=sum(theta)/size(theta)
allocate(xf(lpoints),yf(lpoints),zf(lpoints))
xf(:)=r(:)
!yf(:)=r(:)*theta(:)
!zf(:)=r(:)*sin(theta(:))*phi(:)
yf(:)=rmean*theta(:)
zf(:)=rmean*sin(thetamean)*phi(:)


!GET POSITIONS (CARTESIAN) SET UP FOR MAGNETIC COMPUTATIONS.  THESE ARE PRIMED COORDINATES (SOURCE COORDS, I.E. THE SIM GRID)
allocate(xp(lx1,lx2,lx3),yp(lx1,lx2,lx3),zp(lx1,lx2,lx3))
xp(:,:,:)=x%alt(:,:,:)+Re                               !radial distance from Earth's center
!yp(:,:,:)=xp(:,:,:)*x%theta(:,:,:)                      !southward distance (in the direction of the theta spherical coordinate)
!zp(:,:,:)=xp(:,:,:)*sin(x%theta(:,:,:))*x%phi(:,:,:)    !eastward distance
yp(:,:,:)=rmean*x%theta(:,:,:)
!! the integrations are being treated as Cartesian so flatten out the local spherical coordinates into cartesian, as well
zp(:,:,:)=rmean*sin(thetamean)*x%phi(:,:,:)

!print*, myid2,myid3,'--> field point min/max data:  ',minval(xp),maxval(xp),minval(yp),maxval(yp),minval(zp),maxval(zp)

! differential volumes for source coordinates/integrations
allocate(dV(lx1,lx2,lx3))
allocate(dVend(lx1,lx2),Jxend(lx1,lx2),Jyend(lx1,lx2),Jzend(lx1,lx2))
allocate(Rxend(lx1,lx2),Ryend(lx1,lx2),Rzend(lx1,lx2),Rcubedend(lx1,lx2))
allocate(integrandend(lx1,lx2),integrandavgend(lx1-1,max(lx2-1,1)))

allocate(dVtop(lx1,lx3),Jxtop(lx1,lx3),Jytop(lx1,lx3),Jztop(lx1,lx3))
allocate(Rxtop(lx1,lx3),Rytop(lx1,lx3),Rztop(lx1,lx3),Rcubedtop(lx1,lx3))
allocate(integrandtop(lx1,lx3),integrandavgtop(lx1-1,max(lx3-1,1)))


!> note here that dV's are basically the backward diff volumes; later to be referenced as dV(2:end,2:end,2:end) and so on
if (flag2D/=1) then   !3D differential volume
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        ! avg h's to ix1-1/2, ix2-1/2, ix3-1/2 grid locations
        h1avg=1/8._wp*( x%h1(ix1,ix2,ix3) + x%h1(ix1-1,ix2,ix3) + & 
                   x%h1(ix1,ix2-1,ix3) + x%h1(ix1-1,ix2-1,ix3) + &
                   x%h1(ix1,ix2,ix3-1) + x%h1(ix1-1,ix2,ix3-1) + &
                   x%h1(ix1,ix2-1,ix3-1) + x%h1(ix1-1,ix2-1,ix3-1) ) 
        h2avg=1/8._wp*( x%h2(ix1,ix2,ix3) + x%h2(ix1-1,ix2,ix3) + & 
                   x%h2(ix1,ix2-1,ix3) + x%h2(ix1-1,ix2-1,ix3) + &
                   x%h2(ix1,ix2,ix3-1) + x%h2(ix1-1,ix2,ix3-1) + &
                   x%h2(ix1,ix2-1,ix3-1) + x%h2(ix1-1,ix2-1,ix3-1) ) 
        h3avg=1/8._wp*( x%h3(ix1,ix2,ix3) + x%h3(ix1-1,ix2,ix3) + & 
                   x%h3(ix1,ix2-1,ix3) + x%h3(ix1-1,ix2-1,ix3) + &
                   x%h3(ix1,ix2,ix3-1) + x%h3(ix1-1,ix2,ix3-1) + &
                   x%h3(ix1,ix2-1,ix3-1) + x%h3(ix1-1,ix2-1,ix3-1) ) 
        dV(ix1,ix2,ix3)=h1avg*h2avg*h3avg*x%dx1(ix1)*x%dx2(ix2)*x%dx3(ix3)      !note use here of backward diffs
      end do
    end do
  end do
else                  !plane geometry assumption
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        h1avg=1/4._wp*( x%h1(ix1,ix2,ix3) + x%h1(ix1-1,ix2,ix3) + &
                    x%h1(ix1,ix2,ix3-1) + x%h1(ix1-1,ix2,ix3-1) )
        h3avg=1/4._wp*( x%h3(ix1,ix2,ix3) + x%h3(ix1-1,ix2,ix3) + &
                    x%h3(ix1,ix2,ix3-1) + x%h3(ix1-1,ix2,ix3-1) )
        dV(ix1,ix2,ix3)=h1avg*h3avg*x%dx1(ix1)*x%dx3(ix3)
      end do
    end do
  end do
end if


!> get "end" and "top" pieces for the grid so integrals are not missing any differential volumes
!   The halo_end routine will pass my "begin" and "bottom" pieces of dV to my neighbors on the process grid
call halo_end(dV,dVend,dVtop,tag%dV)
!! need to define the differential volume on the edge of this x3-slab in

!print*, mpi_cfg%myid,' Any problems with dV?  ', any(isnan(dV)),any(isnan(dVend)),any(isnan(dVtop))
!print*, mpi_cfg%myid,'  failures of generate positive dV:  ',any(dV(2:lx1,2:lx2,2:lx3)<=0),any(dVend(2:lx1,2:lx2)<=0), &
!                         any(dVtop(2:lx1,2:lx3)<=0)
!print*, mpi_cfg%myid,'  Min/max edge vals:  ',minval(dVend(2:lx1,2:lx2)),maxval(dVend(2:lx1,2:lx2)), &
!                         minval(dVtop(2:lx1,2:lx3)),maxval(dVend(2:lx1,2:lx3))


!print*, myid2,myid3,'--> dV vals.',minval(dV),maxval(dV),minval(dVend),maxval(dVend),minval(dVtop),maxval(dVtop)

!> Compute projections needed to rotate current density components into magnetic coordinates
allocate(proj_e1er(lx1,lx2,lx3),proj_e2er(lx1,lx2,lx3),proj_e3er(lx1,lx2,lx3))
allocate(proj_e1etheta(lx1,lx2,lx3),proj_e2etheta(lx1,lx2,lx3),proj_e3etheta(lx1,lx2,lx3))
allocate(proj_e1ephi(lx1,lx2,lx3),proj_e2ephi(lx1,lx2,lx3),proj_e3ephi(lx1,lx2,lx3))
allocate(alt(lx1,lx2,lx3))
alt(:,:,:)=x%alt
proj_e1er(:,:,:)=sum(x%e1*x%er,4)
!! fourth dimension of unit vectors is the 3 Cartesian components of each vector
proj_e2er(:,:,:)=sum(x%e2*x%er,4)
proj_e3er(:,:,:)=sum(x%e3*x%er,4)
proj_e1etheta(:,:,:)=sum(x%e1*x%etheta,4)
proj_e2etheta(:,:,:)=sum(x%e2*x%etheta,4)
proj_e3etheta(:,:,:)=sum(x%e3*x%etheta,4)
proj_e1ephi(:,:,:)=sum(x%e1*x%ephi,4)
proj_e2ephi(:,:,:)=sum(x%e2*x%ephi,4)
proj_e3ephi(:,:,:)=sum(x%e3*x%ephi,4)


!> DEALLOCATE GRID MODULE VARIABLES TO SAVE MEMORY (PROGRAM DOESN'T ACTUALLY NEED THESE ONCE X,Y,Z CREATED
call clear_grid(x)
deallocate(r,theta,phi)


!> STORAGE FOR MAGNETIC FIELD CALCULATIONS
allocate(Rx(lx1,lx2,lx3),Ry(lx1,lx2,lx3),Rz(lx1,lx2,lx3))
allocate(Rmag(lx1,lx2,lx3))
allocate(Rcubed(lx1,lx2,lx3))
allocate(integrand(lx1,lx2,lx3),integrandavg(lx1-1,max(lx2-1,1),lx3-1))
!! latter is cell centered hence -1 in size, max is needed to prevent zero sized array
allocate(Br(lpoints),Btheta(lpoints),Bphi(lpoints))
allocate(Brall(lpoints),Bthetaall(lpoints),Bphiall(lpoints))
!! only used by root, but I think workers need to have space allocated for this

! warn user if they are doing a calculation sans parallel current, it's debatable whether this should just error out...
if (mpi_cfg%myid==0 .and. .not. (cfg%flagJpar) ) then
  print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*, 'magcalc WARNING:  you appear to be computing magnetic fields for a simulation with', &
             'parallel currents turned off; results are likely to have substantial error!'
  print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end if


!! MAIN LOOP
!it=1
t=0
tout=t
call dateinc(cfg%dtout,ymd,UTsec)
!! skip first file
it=3
!! don't trigger any special adaptations to filename
main : do while (t < cfg%tdur)
  !TIME STEP CALCULATION
  dt=cfg%dtout    !only compute magnetic field at times when we've done output


  !READ IN THE FULL PLASMA AND FIELD DATA FROM THE OUTPUT FILE (NOTE THAT WE NEED TO KNOW OUTPUT TYPE DONE)
  call input_plasma_currents(cfg%outdir, cfg%out_format, cfg%flagoutput,ymd,UTsec,J1,J2,J3)    !now everyone has their piece of data


  !FORCE PARALLEL CURRENTS TO ZERO BELOW SOME ALTITUDE LIMIT
  if(mpi_cfg%myid==0) print *, 'Zeroing out low altitude currents (these are basically always artifacts)...'
  where (alt < 75000)
    J1=0
    J2=0
    J3=0
  end where


  !DEAL WITH THE WEIRD EDGE ARTIFACTS THAT WE GET IN THE PARALLEL CURRENT
  !SOMETIMES, THIS IS FOR THE X2 DIRECTION
  if(mpi_cfg%myid==0) print *, 'Fixing current edge artifacts...'
  if (mpi_cfg%myid3==mpi_cfg%lid3-1) then
    if (lx3>2) then    !do a ZOH
      J1(:,:,lx3-1)=J1(:,:,lx3-2)
      J1(:,:,lx3)=J1(:,:,lx3-2)
    else
      J1(:,:,lx3-1)=0
      J1(:,:,lx3)=0
    end if
  end if
  if (mpi_cfg%myid3==0) then
    if (lx3>2) then    !do a ZOH
      J1(:,:,1)=J1(:,:,3)
      J1(:,:,2)=J1(:,:,3)
    else
      J1(:,:,1)=0
      J1(:,:,2)=0
    end if
  end if

  !X3 EDGES
  if (mpi_cfg%myid2==mpi_cfg%lid2-1) then
    if (lx2>2) then    !do a ZOH
      J1(:,lx2-1,:)=J1(:,lx2-2,:)
      J1(:,lx2,:)=J1(:,lx2-2,:)
    else
      J1(:,lx2-1,:)=0
      J1(:,lx2-1,:)=0
    end if
  end if
  if (mpi_cfg%myid2==0) then
    if (lx2>2) then    !do a ZOH
      J1(:,1,:)=J1(:,3,:)
      J1(:,2,:)=J1(:,3,:)
    else
      J1(:,1,:)=0
      J1(:,2,:)=0
    end if
  end if


  !ROTATE MAGNETIC FIELDS INTO VERTICAL,SOUTH,EAST COMPONENTS
  if (mpi_cfg%myid==0) then
    print *, 'magcalc.f90 --> Rotating currents into geomagnetic coordinates...'
  end if
  Jx=J1*proj_e1er+J2*proj_e2er+J3*proj_e3er                 !vertical
  Jy=J1*proj_e1etheta+J2*proj_e2etheta+J3*proj_e3etheta     !south
  Jz=J1*proj_e1ephi+J2*proj_e2ephi+J3*proj_e3ephi           !east

!    print *, myid2,myid3,'  --> Min/max values of current',minval(Jx),maxval(Jx),minval(Jy),maxval(Jy), &
!                                               minval(Jz),maxval(Jz)


  !GATHER THE END DATA SO WE DON'T LEAVE OUT A POINT IN THE INTEGRATION
  call halo_end(Jx,Jxend,Jxtop,tag%Jx)
  call halo_end(Jy,Jyend,Jytop,tag%Jy)
  call halo_end(Jz,Jzend,Jztop,tag%Jz)

    !print *, myid2,myid3,'  --> Min/max values of end current',minval(Jxend),maxval(Jxend),minval(Jyend),maxval(Jyend), &
    !                                           minval(Jzend),maxval(Jzend)
    !print *, myid2,myid3,'  --> Min/max values of top current',minval(Jxtop),maxval(Jxtop),minval(Jytop),maxval(Jytop), &
    !                                           minval(Jztop),maxval(Jztop)

  !COMPUTE MAGNETIC FIELDS
  do ipoints=1,lpoints
    if (mpi_cfg%myid == 0 .and. mod(ipoints,100)==0 .and. debug) then
      print *, 'magcalc.f90 --> Computing magnetic field for field point:  ',ipoints,' out of:  ',lpoints
      print *, '            --> ...for time:  ',ymd,UTsec
    end if

    ! compute r-r'
    Rx(:,:,:)=xf(ipoints)-xp(:,:,:)
    Ry(:,:,:)=yf(ipoints)-yp(:,:,:)
    Rz(:,:,:)=zf(ipoints)-zp(:,:,:)
    call halo_end(Rx,Rxend,Rxtop,tag%Rx)
    call halo_end(Ry,Ryend,Rytop,tag%Ry)
    call halo_end(Rz,Rzend,Rztop,tag%Rz)

    ! separately compute average distance for the denominator help with regulation issue and accounts for averaging over each differential volumes
    ! FIXME:  seg fault due to i-1 where i==1
    Rmag=0._wp
    Rmag=1e9     ! will kill off the edges...
    do ix3=2,lx3
      do ix2=2,lx2
        do ix1=2,lx1
          Rmag(ix1,ix2,ix3)=1/8._wp*( sqrt(Rx(ix1,ix2,ix3)**2+Ry(ix1,ix2,ix3)**2+Rz(ix1,ix2,ix3)**2) + &
                                      sqrt(Rx(ix1,ix2,ix3-1)**2+Ry(ix1,ix2,ix3-1)**2+Rz(ix1,ix2,ix3-1)**2) + &
                                      sqrt(Rx(ix1,ix2-1,ix3)**2+Ry(ix1,ix2-1,ix3)**2+Rz(ix1,ix2-1,ix3)**2) + &
                                      sqrt(Rx(ix1,ix2-1,ix3-1)**2+Ry(ix1,ix2-1,ix3-1)**2+Rz(ix1,ix2-1,ix3-1)**2) + &
                                      sqrt(Rx(ix1-1,ix2,ix3)**2+Ry(ix1-1,ix2,ix3)**2+Rz(ix1-1,ix2,ix3)**2) + &
                                      sqrt(Rx(ix1-1,ix2,ix3-1)**2+Ry(ix1-1,ix2,ix3-1)**2+Rz(ix1-1,ix2,ix3-1)**2) + &
                                      sqrt(Rx(ix1-1,ix2-1,ix3)**2+Ry(ix1-1,ix2-1,ix3)**2+Rz(ix1-1,ix2-1,ix3)**2) + &
                                      sqrt(Rx(ix1-1,ix2-1,ix3-1)**2+Ry(ix1-1,ix2-1,ix3-1)**2+Rz(ix1-1,ix2-1,ix3-1)**2) )
        end do
      end do
    end do


    !FIXME:  also need to be averaged using "end" and "top" values for the max x2,x3 edges of this slab...

    !print*, mpi_cfg%myid,' Any problems with Rmag?  ',any(isnan(Rmag))
    !print*, mpi_cfg%myid, 'Min/max vals. for Rmag:  ',minval(Rmag(2:lx1,2:lx2,2:lx3)),maxval(Rmag(2:lx1,2:lx2,2:lx3))


    if (flag2D/=1) then
!      Rcubed(:,:,:)=(Rx**2 + Ry**2 + Rz**2)**(3._wp/2)   !this really is R**3, this has issues with aliasing and div by zero
!      Rcubed(:,:,:)=(Rx**2 + Ry**2 + Rz**2 + Rmin**2)**(3._wp/2)   !this really is R**3, this addresses aliasing but causes underestimate in the value of the computed magnetic fields.  
      Rcubed=Rmag**3

      call halo_end(Rcubed,Rcubedend,Rcubedtop,tag%Rcubed)
      if(mpi_cfg%myid3==mpi_cfg%lid3-1) Rcubedend=R3min     !< avoids div by zero on the end which is set by the haloing
      if(mpi_cfg%myid2==mpi_cfg%lid2-1) Rcubedtop=R3min


      !! FIXME: MAY BE MISSING A CORNER POINT HERE???  NO I THINK IT'S OKAY BASED ON SOME SQUARES I DREW, haha...
      !Bx calculation
      integrand(:,:,:)=mu0/4/pi*(Jy*Rz-Jz*Ry)        ! numerator (to be averaged separately from the denominator
      integrandend(:,:)=mu0/4/pi*(Jyend*Rzend-Jzend*Ryend)
      integrandtop(:,:)=mu0/4/pi*(Jytop*Rztop-Jztop*Rytop)
      integrandavg(:,:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )/Rcubed(2:lx1,2:lx2,2:lx3)

      integrandavgend(:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,lx3) + integrand(2:lx1,1:lx2-1,lx3) + &
                           integrand(1:lx1-1,2:lx2,lx3) + integrand(2:lx1,2:lx2,lx3) )/Rcubed(2:lx1,2:lx2,lx3) + &
                           1/8._wp* (integrandend(1:lx1-1,1:lx2-1) + integrandend(2:lx1,1:lx2-1) + &
                           integrandend(1:lx1-1,2:lx2) + integrandend(2:lx1,2:lx2) )/Rcubedend(2:lx1,2:lx2)

      integrandavgtop(:,:)=1/8._wp*( integrand(1:lx1-1,lx2,1:lx3-1) + integrand(2:lx1,lx2,1:lx3-1) + &
                           integrand(1:lx1-1,lx2,2:lx3) + integrand(2:lx1,lx2,2:lx3) )/Rcubed(2:lx1,lx2,2:lx3) + &
                           1/8._wp*(integrandtop(1:lx1-1,1:lx3-1) + integrandtop(2:lx1,1:lx3-1) + &
                           integrandtop(1:lx1-1,2:lx3) + integrandtop(2:lx1,2:lx3) )/Rcubedtop(2:lx1,2:lx3)

      Br(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))+sum(integrandavgend*dVend(2:lx1,2:lx2))+ &
                    sum(integrandavgtop*dVtop(2:lx1,2:lx3))

!      print*, mpi_cfg%myid,' Br issues:  ',Br(ipoints),any(isnan(Rcubed)),any(isnan(Rcubedend)),any(isnan(Rcubedtop)), &
!               any(isnan(integrandavg)),any(isnan(integrandavgend)),any(isnan(integrandavgtop)),isnan(Br(ipoints))
!      print*, mpi_cfg%myid,' Sum components:  ',isnan(sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))), &
!                                               isnan(sum(integrandavgend*dVend(2:lx1,2:lx2))), &
!                                               isnan(sum(integrandavgtop*dVtop(2:lx1,2:lx3))), any(isnan(dV(2:lx1,2:lx2,2:lx3))), &
!                                               any(isnan(dVend(2:lx1,2:lx2))),any(isnan(dVtop(2:lx1,2:lx3)))
!      print*, mpi_cfg%myid,' Min/max integrandavgs:  ',minval(integrandavg),maxval(integrandavg), &
!                                                       minval(integrandend),maxval(integrandavgend), &
!                                                       minval(integrandavgtop),maxval(integrandavgtop)
!
      !By
      integrand(:,:,:)=-mu0/4/pi*(Jx*Rz-Jz*Rx)
      integrandend(:,:)=-mu0/4/pi*(Jxend*Rzend-Jzend*Rxend)
      integrandtop(:,:)=-mu0/4/pi*(Jxtop*Rztop-Jztop*Rxtop)
      integrandavg(:,:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )/Rcubed(2:lx1,2:lx2,2:lx3)

      integrandavgend(:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,lx3) + integrand(2:lx1,1:lx2-1,lx3) + &
                           integrand(1:lx1-1,2:lx2,lx3) + integrand(2:lx1,2:lx2,lx3) )/Rcubed(2:lx1,2:lx2,lx3) + &
                           1/8._wp*(integrandend(1:lx1-1,1:lx2-1) + integrandend(2:lx1,1:lx2-1) + &
                           integrandend(1:lx1-1,2:lx2) + integrandend(2:lx1,2:lx2) )/Rcubedend(2:lx1,2:lx2)

      integrandavgtop(:,:)=1/8._wp*( integrand(1:lx1-1,lx2,1:lx3-1) + integrand(2:lx1,lx2,1:lx3-1) + &
                           integrand(1:lx1-1,lx2,2:lx3) + integrand(2:lx1,lx2,2:lx3) )/Rcubed(2:lx1,lx2,2:lx3) + &
                           1/8._wp*(integrandtop(1:lx1-1,1:lx3-1) + integrandtop(2:lx1,1:lx3-1) + &
                           integrandtop(1:lx1-1,2:lx3) + integrandtop(2:lx1,2:lx3) )/Rcubedtop(2:lx1,2:lx3)

      Btheta(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))+sum(integrandavgend*dVend(2:lx1,2:lx2))+ &
                        sum(integrandavgtop*dVtop(2:lx1,2:lx3))


      !Bz
      integrand(:,:,:)=mu0/4/pi*(Jx*Ry-Jy*Rx)/Rcubed
      integrandend(:,:)=mu0/4/pi*(Jxend*Ryend-Jyend*Rxend)/Rcubedend
      integrandtop(:,:)=mu0/4/pi*(Jxtop*Rytop-Jytop*Rxtop)/Rcubedtop
      integrandavg(:,:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )/Rcubed(2:lx1,2:lx2,2:lx3)

      integrandavgend(:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,lx3) + integrand(2:lx1,1:lx2-1,lx3) + &
                           integrand(1:lx1-1,2:lx2,lx3) + integrand(2:lx1,2:lx2,lx3) )/Rcubed(2:lx1,2:lx2,lx3) + &
                           1/8._wp*(integrandend(1:lx1-1,1:lx2-1) + integrandend(2:lx1,1:lx2-1) + &
                           integrandend(1:lx1-1,2:lx2) + integrandend(2:lx1,2:lx2) )/Rcubedend(2:lx1,2:lx2)

      integrandavgtop(:,:)=1/8._wp*( integrand(1:lx1-1,lx2,1:lx3-1) + integrand(2:lx1,lx2,1:lx3-1) + &
                           integrand(1:lx1-1,lx2,2:lx3) + integrand(2:lx1,lx2,2:lx3) )/Rcubed(2:lx1,lx2,2:lx3) + &
                           1/8._wp*(integrandtop(1:lx1-1,1:lx3-1) + integrandtop(2:lx1,1:lx3-1) + &
                           integrandtop(1:lx1-1,2:lx3) + integrandtop(2:lx1,2:lx3) )/Rcubedend(2:lx1,2:lx3)

      Bphi(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))+sum(integrandavgend*dVend(2:lx1,2:lx2))+ &
                        sum(integrandavgtop*dVtop(2:lx1,2:lx3))
    else
      Rcubed(:,:,:)=Rx**2+Ry**2    !not really R**3 in 2D, just the denominator of the integrand

!FIXME:  also need a regulator here...
!      where(Rcubed<R3min)
!        Rcubed=R3min    !should be R**2???
!      end where
      call halo_end(Rcubed,Rcubedend,Rcubedtop,tag%Rcubed)
      !! DO WE NEED TO CHECK HERE FOR DIV BY ZERO???
      !! ALSO IN 2D WE KNOW THAT WE ARE ONLY DIVIDED IN THE 3 DIMENSION SO THERE IS NO NEED TO WORRY ABOUT ADDING A 'TOP' ETC.

      !Bx
      integrand(:,:,:)=mu0/4/pi*(-2*Jz*Ry)
      integrandend(:,:)=mu0/4/pi*(-2*Jzend*Ryend)

      integrandavg(:,:,:)=1/4._wp*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx3) + integrand(2:lx1,:,2:lx3) )/Rcubed(2:lx1,:,2:lx3)

      integrandavgend(:,:)=1/4._wp*( integrand(1:lx1-1,:,lx3) + integrand(2:lx1,:,lx3) + &
                           integrandend(1:lx1-1,:) + integrandend(2:lx1,:) )/Rcubedend(2:lx1,:)

      Br(ipoints)=sum(integrandavg*dV(2:lx1,:,2:lx3))+sum(integrandavgend*dVend(2:lx1,:))

      !By
      integrand(:,:,:)=mu0/4/pi*(2*Jz*Rx)
      integrandend(:,:)=mu0/4/pi*(2*Jzend*Rxend)

      integrandavg(:,:,:)=1/4._wp*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx3) + integrand(2:lx1,:,2:lx3) )/Rcubed(2:lx1,:,2:lx3)

      integrandavgend(:,:)=1/4._wp*( integrand(1:lx1-1,:,lx3) + integrand(2:lx1,:,lx3) + &
                           integrandend(1:lx1-1,:) + integrandend(2:lx1,:) )/Rcubedend(2:lx1,:)

      Btheta(ipoints) = sum(integrandavg*dV(2:lx1,:,2:lx3))+sum(integrandavgend*dVend(2:lx1,:))
      !! without dim= input it just sums everything which is what we want

      !Bz
      integrand(:,:,:)=mu0/4/pi*2*(Jx*Ry-Jy*Rx)
      integrandend(:,:)=mu0/4/pi*2*(Jxend*Ryend-Jyend*Rxend)/Rcubedend

      integrandavg(:,:,:)=1/4._wp*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx3) + integrand(2:lx1,:,2:lx3) )/Rcubed(2:lx1,:,2:lx3)

      integrandavgend(:,:)=1/4._wp*( integrand(1:lx1-1,:,lx3) + integrand(2:lx1,:,lx3) + &
                           integrandend(1:lx1-1,:) + integrandend(2:lx1,:) )/Rcubedend(2:lx1,:)

      Bphi(ipoints) = sum(integrandavg*dV(2:lx1,:,2:lx3))+sum(integrandavgend*dVend(2:lx1,:))
      !! without dim= input it just sums everything which is what we want
    end if
  end do


  !A REDUCE OPERATION IS NEEDED HERE TO COMBINE MAGNETIC FIELDS (LINEAR SUPERPOSITION) FROM ALL WORKERS
  if (mpi_cfg%myid ==0) then
    if(debug) print *, 'Attempting reduction of magnetic field...'
  end if
  call mpi_reduce(Br,Brall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(Btheta,Bthetaall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(Bphi,Bphiall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (mpi_cfg%myid == 0) then
    if(debug) print *, 'magcalc.f90 --> Reduced magnetic field...'
    if(debug) print *, '  --> Min/max values of reduced field',minval(Brall),maxval(Brall),minval(Bthetaall),maxval(Bthetaall), &
                                               minval(Bphiall),maxval(Bphiall)
  end if


  if (cfg%dryrun) then
    ierr = mpibreakdown()
    if (ierr /= 0) error stop 'MAGCALC: dry run MPI shutdown failure'
    stop "OK: MAGCALC dry run"
  endif


  !OUTPUT SHOULD BE DONE FOR EVERY INPUT FILE THAT HAS BEEN READ IN
  if (mpi_cfg%myid==0) then
    call cpu_time(tstart)
    call output_magfields(cfg%outdir,ymd,UTsec,Brall,Bthetaall,Bphiall,cfg%out_format)   !mag field data already reduced so just root needs to output
    call cpu_time(tfin)
   if(debug) print *, 'magcalc.f90 --> Output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
  end if
  tout = tout + cfg%dtout


  !NOW OUR SOLUTION IS FULLY UPDATED SO UPDATE TIME VARIABLES TO MATCH...
  it = it + 1
  t = t + dt
  if (mpi_cfg%myid==0 .and. debug) print *, 'magcalc: Moving on to time step (in sec):  ',t,'; end time of simulation:  ',cfg%tdur
  call dateinc(dt,ymd,UTsec)
  if (mpi_cfg%myid==0) print *, 'magcalc.f90 --> Current date',ymd,'Current UT time:  ',UTsec
end do main


!! DEALLOCATE MAIN PROGRAM DATA
deallocate(J1,J2,J3)
deallocate(Jx,Jy,Jz)
deallocate(xp,yp,zp)
deallocate(xf,yf,zf)
deallocate(dV)
deallocate(Rx,Ry,Rz)
deallocate(proj_e1er,proj_e2er,proj_e3er)
deallocate(proj_e1etheta,proj_e2etheta,proj_e3etheta)
deallocate(proj_e1ephi,proj_e2ephi,proj_e3ephi)
deallocate(Rcubed)
deallocate(Rmag)
deallocate(Rcubedend)
deallocate(Rcubedtop)
deallocate(integrand,integrandavg)
deallocate(Br,Btheta,Bphi)
deallocate(dVend,Jxend,Jyend,Jzend,Rxend,Ryend,Rzend)
deallocate(integrandend,integrandavgend)
deallocate(dVtop,Jxtop,Jytop,Jztop,Rxtop,Rytop,Rztop)
deallocate(integrandtop,integrandavgtop)

!! SHUT DOWN MPI
ierr = mpibreakdown()

if (ierr /= 0) then
  write(stderr, *) 'MAGCALC: abnormal MPI shutdown code', ierr, 'Process #', mpi_cfg%myid,' /',mpi_cfg%lid-1
  error stop
endif

print '(A)', 'MAGCALC: complete'

end program
