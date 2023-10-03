Program MagCalc
!! Need program statement for FORD
!!
!! This program computes magnetic field fluctuations using current density from a gemini output simulation

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use magcalc_cli, only : cli
use phys_consts, only : pi,mu0, wp, re, debug
use grid, only : lx1, lx2, lx3, lx2all,lx3all,grid_size,read_grid,calc_subgrid_size
use meshobj, only : curvmesh
use timeutils, only : dateinc,find_time_elapsed
use gemini3d_config, only : gemini_cfg
use io, only : input_plasma_currents,create_outdir_mag,output_magfields
use mpimod, only: mpibreakdown, process_grid_auto, mpi_manualgrid, halo_end, &
mpi_cfg, mpi_realprec, tag=>gemini_mpi
use h5fortran, only : hdf5_file
use filesystem, only : suffix

use mpi_f08, only: mpi_init,mpi_finalize,mpi_comm_rank,mpi_reduce,mpi_sum, mpi_comm_world

implicit none (type, external)

!> VARIABLES READ IN FROM CONFIG FILE

real(wp) :: UTsec
!! UT (s)
integer, dimension(3) :: ymd
!! year, month, day (current)

type(gemini_cfg) :: cfg
!! most user parameters

!! GRID STRUCTURE
class(curvmesh), pointer :: x
!! structure containing grid locations, finite differences, etc.:  see grid module for details

!STATE VARIABLES
real(wp), dimension(:,:,:), allocatable :: J1,J2,J3      !electrodynamic state variables

!TEMPORAL VARIABLES
real(wp) :: t=0, dt      !time from beginning of simulation (s) and time step (s)
real(wp) :: tout    !time for next output and time between outputs
real(wp) :: tstart,tfin   !temp. vars. for measuring performance of code blocks
integer :: it        !time and species loop indices

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
real(wp), dimension(:,:,:), allocatable :: integrand
real(wp), dimension(:,:,:), allocatable :: alt
real(wp), dimension(:), allocatable :: Br,Btheta,Bphi
real(wp), dimension(:), allocatable :: Brall,Bthetaall,Bphiall
real(wp), dimension(:,:), allocatable :: Jxend,Jyend,Jzend,Rxend,Ryend,Rzend,Rcubedend,dVend,Rmagend
real(wp), dimension(:,:), allocatable :: integrandend

real(wp), dimension(:,:), allocatable :: Jxtop,Jytop,Jztop,Rxtop,Rytop,Rztop,Rcubedtop,dVtop,Rmagtop
real(wp), dimension(:,:), allocatable :: integrandtop
real(wp), dimension(:), allocatable :: integrandcorner

real(wp), dimension(:,:), allocatable :: xpend,ypend,zpend
real(wp), dimension(:,:), allocatable :: xptop,yptop,zptop

real(wp), dimension(:), allocatable :: dVcorner,xpcorner,ypcorner,zpcorner,Jxcorner,Jycorner,Jzcorner,Rmagcorner,Rcubedcorner
real(wp), dimension(:), allocatable :: Rxcorner,Rycorner,Rzcorner

integer :: ix1,ix2,ix3
real(wp) :: rmean,thetamean

!! FOR SPECIFYING THE PROCESS GRID
integer :: lid2in,lid3in

!! REGULATOR FOR 1/R^3
!real(wp), parameter :: R3min=1d11     !works well for 192x192 in situ
!real(wp), parameter :: R3min=1d9
!real(wp), parameter :: Rmin=5d3

!! for keeping track of start and end times requested by the user
integer, dimension(3) :: ymdstart,ymdend,ymdfinal
real(wp) :: UTsecstart,UTsecend,telend,UTsecfinal
real(wp) :: h1avg,h2avg,h3avg
real(wp), dimension(:,:,:), allocatable :: Rmag

!! --- MAIN PROGRAM
call mpi_init()

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
call calc_subgrid_size(lx2all,lx3all)
print*, 'grid size:  ',lx1,lx2,lx3,lx2all,lx3all

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
allocate(J1(-1:lx1+2,-1:lx2+2,-1:lx3+2),J2(-1:lx1+2,-1:lx2+2,-1:lx3+2),J3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(Jx(lx1,lx2,lx3),Jy(lx1,lx2,lx3),Jz(lx1,lx2,lx3))


!NOW DEAL WITH THE UNPRIMED COORDINATES
block
  integer :: u
  type(hdf5_file) :: hf

  select case (suffix(cfg%indatsize))
  case ('.dat')
    open(newunit=u,file=cfg%fieldpointfile,status='old',form='unformatted',access='stream',action='read')
    read(u) lpoints    !size of coordinates for field points
    if (mpi_cfg%myid==0) print *, 'magcalc.f90 --> Number of field points:  ',lpoints
    allocate(r(lpoints),theta(lpoints),phi(lpoints))
    read(u) r,theta,phi
    close(u)
  case ('.h5')
    !! hdf5 file input
    call hf%open(cfg%fieldpointfile, action='r')

    call hf%read('/lpoints',lpoints)
    allocate(r(lpoints), theta(lpoints), phi(lpoints))
    call hf%read('/r',r)
    call hf%read('/theta',theta)
    call hf%read('/phi',phi)
    call hf%close()
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
if (mpi_cfg%myid==0) print*, 'magcalc.f90 --> setting up field point x,y,z...'
allocate(xp(lx1,lx2,lx3),yp(lx1,lx2,lx3),zp(lx1,lx2,lx3))
xp(:,:,:)=x%alt(1:lx1,1:lx2,1:lx3)+Re                               !radial distance from Earth's center
!yp(:,:,:)=xp(:,:,:)*x%theta(:,:,:)                      !southward distance (in the direction of the theta spherical coordinate)
!zp(:,:,:)=xp(:,:,:)*sin(x%theta(:,:,:))*x%phi(:,:,:)    !eastward distance
yp(:,:,:)=rmean*x%theta(1:lx1,1:lx2,1:lx3)
!! the integrations are being treated as Cartesian so flatten out the local spherical coordinates into cartesian, as well
zp(:,:,:)=rmean*sin(thetamean)*x%phi(1:lx1,1:lx2,1:lx3)

!print*, myid2,myid3,'--> field point min/max data:  ',minval(xp),maxval(xp),minval(yp),maxval(yp),minval(zp),maxval(zp)

! differential volumes for source coordinates/integrations
if (mpi_cfg%myid==0) print*, 'magcalc.f90 --> computing differential volumes for integral(s)...'
allocate(dV(lx1,lx2,lx3))
allocate(dVend(lx1,lx2),Jxend(lx1,lx2),Jyend(lx1,lx2),Jzend(lx1,lx2))
allocate(Rxend(lx1,lx2),Ryend(lx1,lx2),Rzend(lx1,lx2),Rcubedend(lx1,lx2),Rmagend(lx1,lx2))
allocate(integrandend(lx1,lx2))

allocate(dVtop(lx1,lx3),Jxtop(lx1,lx3),Jytop(lx1,lx3),Jztop(lx1,lx3))
allocate(Rxtop(lx1,lx3),Rytop(lx1,lx3),Rztop(lx1,lx3),Rcubedtop(lx1,lx3),Rmagtop(lx1,lx3))
allocate(integrandtop(lx1,lx3))

allocate(xpend(lx1,lx2),ypend(lx1,lx2),zpend(lx1,lx2))
allocate(xptop(lx1,lx3),yptop(lx1,lx3),zptop(lx1,lx3))

allocate(dVcorner(lx1),xpcorner(lx1),ypcorner(lx1),zpcorner(lx1),Jxcorner(lx1),Jycorner(lx1),Jzcorner(lx1), &
           Rmagcorner(lx1),Rcubedcorner(lx1))
allocate(Rxcorner(lx1),Rycorner(lx1),Rzcorner(lx1))
allocate(integrandcorner(lx1))

!> note here that dV's are basically the backward diff volumes; later to be referenced as dV(2:end,2:end,2:end) and so on.
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


if (mpi_cfg%myid==0) print*, 'magcalc.f90 --> worker exchange of edge volumes...'
! FIXME: does this need message passing???  Seems like these coudl be computed locally since the ghost cell metric factors and differentials are already stored by workers...
!> get "end" and "top" pieces for the grid so integrals are not missing any differential volumes
!   The halo_end routine will pass my "begin" and "bottom" pieces of dV to neighbors on the process grid
call halo_end(dV,dVend,dVtop,dVcorner,tag%dV)
!! need to define the differential volume on the edge of this x3-slab in

if (mpi_cfg%myid==0) print*, 'magcalc.f90 --> worker exchange of edge distances...'
!> now get the "end" and "top" pieces for the source coordinates
call halo_end(xp,xpend,xptop,xpcorner,tag%Rx)    !just reuse position tag
call halo_end(yp,ypend,yptop,ypcorner,tag%Ry)
call halo_end(zp,zpend,zptop,zpcorner,tag%Rz)

!> Compute projections needed to rotate current density components into magnetic coordinates
if (mpi_cfg%myid==0) print*, 'magcalc.f90 --> calculating projections need for rotation of current densities...'
allocate(proj_e1er(lx1,lx2,lx3),proj_e2er(lx1,lx2,lx3),proj_e3er(lx1,lx2,lx3))
allocate(proj_e1etheta(lx1,lx2,lx3),proj_e2etheta(lx1,lx2,lx3),proj_e3etheta(lx1,lx2,lx3))
allocate(proj_e1ephi(lx1,lx2,lx3),proj_e2ephi(lx1,lx2,lx3),proj_e3ephi(lx1,lx2,lx3))
allocate(alt(lx1,lx2,lx3))
alt(1:lx1,1:lx2,1:lx3)=x%alt
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
!call clear_grid(x)
deallocate(r,theta,phi)


!> STORAGE FOR MAGNETIC FIELD CALCULATIONS
allocate(Rx(lx1,lx2,lx3),Ry(lx1,lx2,lx3),Rz(lx1,lx2,lx3))
allocate(Rmag(lx1,lx2,lx3))
allocate(Rcubed(lx1,lx2,lx3))
allocate(integrand(lx1,lx2,lx3))
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
  if (mpi_cfg%myid==0) then
    print *, 'magcalc.f90 --> Reading input datafile for time:  ',ymd,UTsec
  end if
  call input_plasma_currents(cfg%outdir, cfg%out_format, cfg%flagoutput,ymd,UTsec,J1,J2,J3)    !now everyone has their piece of data

  !! FAC can often have edge artifacts due to boundary being too close to the disturbance being modeled.
  call fixJ(J1,J2,J3)

  !ROTATE MAGNETIC FIELDS INTO VERTICAL,SOUTH,EAST COMPONENTS
  if (mpi_cfg%myid==0) then
    print *, 'magcalc.f90 --> Rotating currents into geomagnetic coordinates...'
  end if
  Jx=J1(1:lx1,1:lx2,1:lx3)*proj_e1er+J2(1:lx1,1:lx2,1:lx3)*proj_e2er+J3(1:lx1,1:lx2,1:lx3)*proj_e3er                 !vertical
  Jy=J1(1:lx1,1:lx2,1:lx3)*proj_e1etheta+J2(1:lx1,1:lx2,1:lx3)*proj_e2etheta+J3(1:lx1,1:lx2,1:lx3)*proj_e3etheta     !south
  Jz=J1(1:lx1,1:lx2,1:lx3)*proj_e1ephi+J2(1:lx1,1:lx2,1:lx3)*proj_e2ephi+J3(1:lx1,1:lx2,1:lx3)*proj_e3ephi           !east

!    print *, myid2,myid3,'  --> Min/max values of current',minval(Jx),maxval(Jx),minval(Jy),maxval(Jy), &
!                                               minval(Jz),maxval(Jz)
  !GATHER THE END DATA SO WE DON'T LEAVE OUT A POINT IN THE INTEGRATION
  call halo_end(Jx,Jxend,Jxtop,Jxcorner,tag%Jx)
  call halo_end(Jy,Jyend,Jytop,Jycorner,tag%Jy)
  call halo_end(Jz,Jzend,Jztop,Jzcorner,tag%Jz)

  !COMPUTE MAGNETIC FIELDS
  do ipoints=1,lpoints
    if (mpi_cfg%myid == 0 .and. mod(ipoints,100)==0 .and. debug) then
      print *, 'magcalc.f90 --> Computing magnetic field for field point:  ',ipoints,' out of:  ',lpoints
      print *, '            --> ...for time:  ',ymd,UTsec
    end if

    ! compute r-r', including at endpoints needed to fully cover coordinate calculations for all differential volumes
    ! note that the primed locations can be computed ONCE at the beginning of the simulation which will vastly reduce the amount of message passing
    Rx(:,:,:)=xf(ipoints)-xp(:,:,:)
    Ry(:,:,:)=yf(ipoints)-yp(:,:,:)
    Rz(:,:,:)=zf(ipoints)-zp(:,:,:)
    Rxend(:,:)=xf(ipoints)-xpend(:,:); Rxtop(:,:)=xf(ipoints)-xptop(:,:)
    Ryend(:,:)=yf(ipoints)-ypend(:,:); Rytop(:,:)=yf(ipoints)-yptop(:,:)
    Rzend(:,:)=zf(ipoints)-zpend(:,:); Rztop(:,:)=zf(ipoints)-zptop(:,:)
    Rxcorner(:)=xf(ipoints)-xpcorner(:)
    Rycorner(:)=yf(ipoints)-ypcorner(:)
    Rzcorner(:)=zf(ipoints)-zpcorner(:)
    call calcRmag(Rx,Ry,Rz,Rxend,Ryend,Rzend,Rxtop,Rytop,Rztop,Rxcorner,Rycorner,Rzcorner, &
                    Rmag,Rmagend,Rmagtop,Rmagcorner)

    if (flag2D/=1) then
      Rcubed=Rmag**3
      Rcubedend=Rmagend**3
      Rcubedtop=Rmagtop**3
      Rcubedcorner=Rmagcorner**3


      !! FIXME: MAY BE MISSING A CORNER POINT HERE???  NO I THINK IT'S OKAY BASED ON SOME SQUARES I DREW, haha...
      ! FIXME:  not okay...
      !Bx calculation
      integrand(:,:,:)=mu0/4/pi*(Jy*Rz-Jz*Ry)        ! numerator (to be averaged separately from the denominator
      integrandend(:,:)=mu0/4/pi*(Jyend*Rzend-Jzend*Ryend)
      integrandtop(:,:)=mu0/4/pi*(Jytop*Rztop-Jztop*Rytop)
      integrandcorner(:)=mu0/4/pi*(Jycorner*Rzcorner-Jzcorner*Rycorner)
      Br(ipoints)=integrate3D(integrand,integrandend,integrandtop,integrandcorner)


      !By
      integrand(:,:,:)=-mu0/4/pi*(Jx*Rz-Jz*Rx)
      integrandend(:,:)=-mu0/4/pi*(Jxend*Rzend-Jzend*Rxend)
      integrandtop(:,:)=-mu0/4/pi*(Jxtop*Rztop-Jztop*Rxtop)
      integrandcorner(:)=-mu0/4/pi*(Jxcorner*Rzcorner-Jzcorner*Rxcorner)
      Btheta(ipoints)=integrate3D(integrand,integrandend,integrandtop,integrandcorner)


      !Bz
      integrand(:,:,:)=mu0/4/pi*(Jx*Ry-Jy*Rx)
      integrandend(:,:)=mu0/4/pi*(Jxend*Ryend-Jyend*Rxend)
      integrandtop(:,:)=mu0/4/pi*(Jxtop*Rytop-Jytop*Rxtop)
      integrandcorner(:)=mu0/4/pi*(Jxcorner*Rycorner-Jycorner*Rxcorner)
      Bphi(ipoints)=integrate3D(integrand,integrandend,integrandtop,integrandcorner)
    else
      Rcubed(:,:,:)=Rx**2+Ry**2    !not really R**3 in 2D, just the denominator of the integrand

!FIXME:  also need a regulator here...
!      where(Rcubed<R3min)
!        Rcubed=R3min    !should be R**2???
!      end where
      call halo_end(Rcubed,Rcubedend,Rcubedtop,Rcubedcorner,tag%Rcubed)    ! corner not used here...
      !! DO WE NEED TO CHECK HERE FOR DIV BY ZERO???
      !! ALSO IN 2D WE KNOW THAT WE ARE ONLY DIVIDED IN THE 3 DIMENSION SO THERE IS NO NEED TO WORRY ABOUT ADDING A 'TOP' ETC.

      !Bx
      integrand(:,:,:)=mu0/4/pi*(-2*Jz*Ry)
      integrandend(:,:)=mu0/4/pi*(-2*Jzend*Ryend)
      Br(ipoints)=integrate2D(integrand,integrandend)

      !By
      integrand(:,:,:)=mu0/4/pi*(2*Jz*Rx)
      integrandend(:,:)=mu0/4/pi*(2*Jzend*Rxend)
      Btheta(ipoints)=integrate2D(integrand,integrandend)

      !Bz
      integrand(:,:,:)=mu0/4/pi*2*(Jx*Ry-Jy*Rx)
      integrandend(:,:)=mu0/4/pi*2*(Jxend*Ryend-Jyend*Rxend)/Rcubedend
      Bphi(ipoints) = integrate2D(integrand,integrandend)
    end if
  end do


  !A REDUCE OPERATION IS NEEDED HERE TO COMBINE MAGNETIC FIELDS (LINEAR SUPERPOSITION) FROM ALL WORKERS
  if (mpi_cfg%myid ==0) then
    if(debug) print *, 'Attempting reduction of magnetic field...'
  end if
  call mpi_reduce(Br,Brall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD)
  call mpi_reduce(Btheta,Bthetaall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD)
  call mpi_reduce(Bphi,Bphiall,lpoints,mpi_realprec,MPI_SUM,0,MPI_COMM_WORLD)
  if (mpi_cfg%myid == 0) then
    if(debug) print *, 'magcalc.f90 --> Reduced magnetic field...'
    if(debug) print *, '  --> Min/max values of reduced field',minval(Brall),maxval(Brall),minval(Bthetaall),maxval(Bthetaall), &
                                               minval(Bphiall),maxval(Bphiall)
  end if


  if (cfg%dryrun) then
    if (mpibreakdown() /= 0) error stop 'MAGCALC: dry run MPI shutdown failure'
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
deallocate(integrand)
deallocate(Br,Btheta,Bphi)
deallocate(dVend,Jxend,Jyend,Jzend,Rxend,Ryend,Rzend)
deallocate(integrandend)
deallocate(dVtop,Jxtop,Jytop,Jztop,Rxtop,Rytop,Rztop)
deallocate(integrandtop)
deallocate(dVcorner,xpcorner,ypcorner,zpcorner,Jxcorner,Jycorner,Jzcorner, &
           Rmagcorner,Rcubedcorner)
deallocate(Rxcorner,Rycorner,Rzcorner)
deallocate(integrandcorner)


!! SHUT DOWN MPI

if (mpibreakdown() /= 0) then
  write(stderr, *) 'MAGCALC: abnormal MPI shutdown: Process #', mpi_cfg%myid,' /',mpi_cfg%lid-1
  error stop
endif

print '(A)', 'MAGCALC: complete'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains    ! declare integral functions as internal subprograms; too specific to be used elsewhere.  Also they access data from the main program unti because I don't feel like including these are arguments.
  subroutine fixJ(J1,J2,J3)
    ! host program data used (but not modified)
    ! mpi_cfg, alt
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: J1,J2,J3
    integer :: lx1,lx2,lx3

    lx1=size(J1,1)-4; lx2=size(J2,2)-4; lx3=size(J3,3)-4;

    !FORCE PARALLEL CURRENTS TO ZERO BELOW SOME ALTITUDE LIMIT
    if(mpi_cfg%myid==0) print *, 'Zeroing out low altitude currents (these are basically always artifacts)...'
    where (x%alt < 75000)
      J1=0
      J2=0
      J3=0
    end where

    !! FAC can often have edge artifacts due to boundary being too close to the disturbance being modeled.
    !  this code will remove these.
    ! x3 global grid edges
    if(mpi_cfg%myid==0) print *, 'Fixing current edge artifacts...'
    if (mpi_cfg%myid3==mpi_cfg%lid3-1) then
      if (lx3>2) then    !do a ZOH from the 3rd to last cell
        J1(:,:,lx3-1)=J1(:,:,lx3-2)
        J1(:,:,lx3)=J1(:,:,lx3-2)
      else
        J1(:,:,lx3-1)=0
        J1(:,:,lx3)=0
      end if
    end if
    if (mpi_cfg%myid3==0) then
      if (lx3>2) then    !do a ZOH from the 3rd cell
        J1(:,:,1)=J1(:,:,3)
        J1(:,:,2)=J1(:,:,3)
      else
        J1(:,:,1)=0
        J1(:,:,2)=0
      end if
    end if

    ! x2 global grid edges
    !if(mpi_cfg%myid==0) print *, 'Fixing current edge artifacts, x3...'
    if (mpi_cfg%myid2==mpi_cfg%lid2-1) then
      if (lx2>2) then
        J1(:,lx2-1,:)=J1(:,lx2-2,:)
        J1(:,lx2,:)=J1(:,lx2-2,:)
      else
        J1(:,lx2-1,:)=0
        J1(:,lx2-1,:)=0
      end if
    end if
    if (mpi_cfg%myid2==0) then
      if (lx2>2) then
        J1(:,1,:)=J1(:,3,:)
        J1(:,2,:)=J1(:,3,:)
      else
        J1(:,1,:)=0
        J1(:,2,:)=0
      end if
    end if
  end subroutine fixJ

  subroutine calcRmag(Rx,Ry,Rz, &
                    Rxend,Ryend,Rzend, &
                    Rxtop,Rytop,Rztop, &
                    Rxcorner,Rycorner,Rzcorner, &
                    Rmag,Rmagend,Rmagtop,Rmagcorner)
    ! host program data used (but not modified)
    ! mpi_cfg
    real(wp), dimension(:,:,:), intent(in) :: Rx,Ry,Rz
    real(wp), dimension(:,:), intent(in) :: Rxend,Ryend,Rzend,Rxtop,Rytop,Rztop
    real(wp), dimension(:), intent(in) :: Rxcorner,Rycorner,Rzcorner
    real(wp), dimension(:,:,:),intent(inout) :: Rmag
    !! intent(out)
    real(wp), dimension(:,:),intent(inout) :: Rmagend,Rmagtop
    !! intent(out)
    real(wp), dimension(:),intent(inout) :: Rmagcorner
    !! intent(out)
    integer :: ix1,ix2,ix3,lx1,lx2,lx3

    lx1=size(Rx,1); lx2=size(Rx,2); lx3=size(Rx,3);

    !! separately compute average distance for the denominator help with regulation issue and
    !! accounts for averaging over each differential volumes
    Rmag = 0
    do ix3=2,lx3
      do ix2=2,lx2
        do ix1=2,lx1
          Rmag(ix1,ix2,ix3)=1/8._wp*( sqrt(Rx(ix1,ix2,ix3)**2      +  Ry(ix1,ix2,ix3)**2      +  Rz(ix1,ix2,ix3)**2) + &       ! i,j,k
                                      sqrt(Rx(ix1,ix2,ix3-1)**2    +  Ry(ix1,ix2,ix3-1)**2    +  Rz(ix1,ix2,ix3-1)**2) + &     ! i,j,k-1
                                      sqrt(Rx(ix1,ix2-1,ix3)**2    +  Ry(ix1,ix2-1,ix3)**2    +  Rz(ix1,ix2-1,ix3)**2) + &     ! i,j-1,k
                                      sqrt(Rx(ix1,ix2-1,ix3-1)**2  +  Ry(ix1,ix2-1,ix3-1)**2  +  Rz(ix1,ix2-1,ix3-1)**2) + &   ! i,j-1,k-1
                                      sqrt(Rx(ix1-1,ix2,ix3)**2    +  Ry(ix1-1,ix2,ix3)**2    +  Rz(ix1-1,ix2,ix3)**2) + &     ! i-1,j,k
                                      sqrt(Rx(ix1-1,ix2,ix3-1)**2  +  Ry(ix1-1,ix2,ix3-1)**2  +  Rz(ix1-1,ix2,ix3-1)**2) + &   ! i-1,j,k-1
                                      sqrt(Rx(ix1-1,ix2-1,ix3)**2  +  Ry(ix1-1,ix2-1,ix3)**2  +  Rz(ix1-1,ix2-1,ix3)**2) + &   ! i-1,j-1,k
                                      sqrt(Rx(ix1-1,ix2-1,ix3-1)**2+  Ry(ix1-1,ix2-1,ix3-1)**2+  Rz(ix1-1,ix2-1,ix3-1)**2) )   ! i-1,j-1,k-1
        end do
      end do
    end do

    ! end and top values should be added.
    Rmagend = 0
    if (mpi_cfg%myid3/=mpi_cfg%lid3-1) then
      do ix2=2,lx2
        do ix1=2,lx1
          Rmagend(ix1,ix2)=1/8._wp*( sqrt(Rx(ix1,ix2,lx3)**2+Ry(ix1,ix2,lx3)**2+Rz(ix1,ix2,lx3)**2) + &
                                        sqrt(Rxend(ix1,ix2)**2+Ryend(ix1,ix2)**2+Rzend(ix1,ix2)**2) + &
                                        sqrt(Rx(ix1,ix2-1,lx3)**2+Ry(ix1,ix2-1,lx3)**2+Rz(ix1,ix2-1,lx3)**2) + &
                                        sqrt(Rxend(ix1,ix2-1)**2+Ryend(ix1,ix2-1)**2+Rzend(ix1,ix2-1)**2) + &
                                        sqrt(Rx(ix1-1,ix2,lx3)**2+Ry(ix1-1,ix2,lx3)**2+Rz(ix1-1,ix2,lx3)**2) + &
                                        sqrt(Rxend(ix1-1,ix2)**2+Ryend(ix1-1,ix2)**2+Rzend(ix1-1,ix2)**2) + &
                                        sqrt(Rx(ix1-1,ix2-1,lx3)**2+Ry(ix1-1,ix2-1,lx3)**2+Rz(ix1-1,ix2-1,lx3)**2) + &
                                        sqrt(Rxend(ix1-1,ix2-1)**2+Ryend(ix1-1,ix2-1)**2+Rzend(ix1-1,ix2-1)**2) )
        end do
      end do
    end if
    Rmagtop = 0
    if (mpi_cfg%myid2/=mpi_cfg%lid2-1) then
      do ix3=2,lx3
        do ix1=2,lx1
          Rmagtop(ix1,ix3)=1/8._wp*( sqrt(Rx(ix1,lx2,ix3)**2+Ry(ix1,lx2,ix3)**2+Rz(ix1,lx2,ix3)**2) + &
                                        sqrt(Rx(ix1,lx2,ix3-1)**2+Ry(ix1,lx2,ix3-1)**2+Rz(ix1,lx2,ix3-1)**2) + &
                                        sqrt(Rxtop(ix1,ix3)**2+Rytop(ix1,ix3)**2+Rztop(ix1,ix3)**2) + &
                                        sqrt(Rxtop(ix1,ix3-1)**2+Rytop(ix1,ix3-1)**2+Rztop(ix1,ix3-1)**2) + &
                                        sqrt(Rx(ix1-1,lx2,ix3)**2+Ry(ix1-1,lx2,ix3)**2+Rz(ix1-1,lx2,ix3)**2) + &
                                        sqrt(Rx(ix1-1,lx2,ix3-1)**2+Ry(ix1-1,lx2,ix3-1)**2+Rz(ix1-1,lx2,ix3-1)**2) + &
                                        sqrt(Rxtop(ix1-1,ix3)**2+Rytop(ix1-1,ix3)**2+Rztop(ix1-1,ix3)**2) + &
                                        sqrt(Rxtop(ix1-1,ix3-1)**2+Rytop(ix1-1,ix3-1)**2+Rztop(ix1-1,ix3-1)**2) )
        end do
      end do
    end if

    ! corner cell distance to be computed
    Rmagcorner = 0
    if (mpi_cfg%myid3/=mpi_cfg%lid3-1 .and. mpi_cfg%myid2/=mpi_cfg%lid2-1) then
      do ix1=2,lx1
        Rmagcorner(ix1)=1/8._wp*( sqrt(Rxcorner(ix1)**2      +  Rycorner(ix1)**2      +  Rzcorner(ix1)**2) + &                    ! i,j,k
                                  sqrt(Rxtop(ix1,lx3)**2    +  Rytop(ix1,lx3)**2    +  Rztop(ix1,lx3)**2) + &                ! i,j,k-1
                                  sqrt(Rxend(ix1,lx2)**2    +  Ryend(ix1,lx2)**2    +  Rzend(ix1,lx2)**2) + &                ! i,j-1,k
                                  sqrt(Rx(ix1,lx2,lx3)**2  +  Ry(ix1,lx2,lx3)**2  +  Rz(ix1,lx2,lx3)**2) + &                 ! i,j-1,k-1
                                  sqrt(Rxcorner(ix1-1)**2    +  Rycorner(ix1-1)**2    +  Rzcorner(ix1-1)**2) + &             ! i-1,j,k
                                  sqrt(Rxtop(ix1-1,lx3)**2  +  Rytop(ix1-1,lx3)**2  +  Rztop(ix1-1,lx3)**2) + &              ! i-1,j,k-1
                                  sqrt(Rxend(ix1-1,lx2)**2  +  Ryend(ix1-1,lx2)**2  +  Rzend(ix1-1,lx2)**2) + &              ! i-1,j-1,k
                                  sqrt(Rx(ix1-1,lx2,lx3)**2+  Ry(ix1-1,lx2,lx3)**2+  Rz(ix1-1,lx2,lx3)**2) )                 ! i-1,j-1,k-1
      end do
    end if
  end subroutine calcRmag


  function integrate3D(integrand,integrandend,integrandtop,integrandcorner)
    ! host program data used (but not modified)
    ! mpi_cfg, Rcubed*, dV*
    real(wp), dimension(:,:,:), intent(in) :: integrand
    real(wp), dimension(:,:), intent(in) :: integrandend
    real(wp), dimension(:,:), intent(in) :: integrandtop
    real(wp), dimension(:), intent(in) :: integrandcorner
    real(wp), dimension(:,:,:), allocatable :: integrandavg
    real(wp), dimension(:,:), allocatable :: integrandavgend,integrandavgtop
    real(wp), dimension(:), allocatable :: integrandavgcorner
    integer :: lx1,lx2,lx3
    real(wp) :: integrate3D

    lx1=size(integrand,1); lx2=size(integrand,2); lx3=size(integrand,3);

    allocate(integrandavg(lx1-1,max(lx2-1,1),max(lx3-1,1)), &
             integrandavgend(lx1-1,max(lx2-1,1)), &
             integrandavgtop(lx1-1,max(lx3-1,1)), &
             integrandavgcorner(lx1-1) )

    integrandavg(:,:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                                  integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                                  integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                                  integrand(1:lx1-1,2:lx2,2:lx3)   + integrand(2:lx1,2:lx2,2:lx3) )/ &
                                 Rcubed(2:lx1,2:lx2,2:lx3)

    integrandavgend = 0
    if (mpi_cfg%myid3/=mpi_cfg%lid3-1) then
      integrandavgend(:,:)=1/8._wp*( integrand(1:lx1-1,1:lx2-1,lx3) + integrand(2:lx1,1:lx2-1,lx3) + &
                           integrand(1:lx1-1,2:lx2,lx3) + integrand(2:lx1,2:lx2,lx3) + &
                           integrandend(1:lx1-1,1:lx2-1) + integrandend(2:lx1,1:lx2-1) + &
                           integrandend(1:lx1-1,2:lx2) + integrandend(2:lx1,2:lx2) )/Rcubedend(2:lx1,2:lx2)
    end if

    integrandavgtop = 0
    if (mpi_cfg%myid2/=mpi_cfg%lid2-1) then
      integrandavgtop(:,:)=1/8._wp*( integrand(1:lx1-1,lx2,1:lx3-1) + integrand(2:lx1,lx2,1:lx3-1) + &
                                     integrand(1:lx1-1,lx2,2:lx3)   + integrand(2:lx1,lx2,2:lx3) + &
                                     integrandtop(1:lx1-1,1:lx3-1) + integrandtop(2:lx1,1:lx3-1) + &
                                     integrandtop(1:lx1-1,2:lx3)   + integrandtop(2:lx1,2:lx3) )/ &
                                    Rcubedtop(2:lx1,2:lx3)
    end if

    integrandavgcorner = 0
    if (mpi_cfg%myid3/=mpi_cfg%lid3-1 .and. mpi_cfg%myid2/=mpi_cfg%lid2-1) then
      integrandavgcorner(:)=1/8._wp*( integrandcorner(1:lx1-1) + integrandcorner(2:lx1) + &
                                     integrandend(1:lx1-1,lx2) + integrandend(2:lx1,lx2) + &
                                     integrandtop(1:lx1-1,lx3) + integrandtop(2:lx1,lx3) + &
                                     integrand(1:lx1-1,lx2,lx3) + integrand(2:lx1,lx2,lx3) &
                                     )/Rcubedcorner(2:lx1)
    end if

    integrate3D=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))+sum(integrandavgend*dVend(2:lx1,2:lx2))+ &
                      sum(integrandavgtop*dVtop(2:lx1,2:lx3))+sum(integrandavgcorner*dVcorner(2:lx1))

    deallocate(integrandavg,integrandavgend,integrandavgtop,integrandavgcorner)
  end function integrate3D

  function integrate2D(integrand,integrandend)
    ! host program data used (but not modified)
    ! mpi_cfg, Rcubed*, dV*
    real(wp), dimension(:,:,:), intent(in) :: integrand
    real(wp), dimension(:,:), intent(in) :: integrandend
    real(wp), dimension(:,:,:), allocatable :: integrandavg
    real(wp), dimension(:,:), allocatable :: integrandavgend
    integer :: lx1,lx2,lx3
    real(wp) :: integrate2D

    lx1=size(integrand,1); lx2=size(integrand,2); lx3=size(integrand,3);

    allocate(integrandavg(lx1-1,max(lx2-1,1),max(lx3-1,1)), &
             integrandavgend(lx1-1,max(lx2-1,1)) )

    integrandavg(:,:,:)=1/4._wp*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                         integrand(1:lx1-1,:,2:lx3) + integrand(2:lx1,:,2:lx3) )/Rcubed(2:lx1,:,2:lx3)

    integrandavgend(:,:)=1/4._wp*( integrand(1:lx1-1,:,lx3) + integrand(2:lx1,:,lx3) + &
                         integrandend(1:lx1-1,:) + integrandend(2:lx1,:) )/Rcubedend(2:lx1,:)

    integrate2D=sum(integrandavg*dV(2:lx1,:,2:lx3))+sum(integrandavgend*dVend(2:lx1,:))

    deallocate(integrandavg,integrandavgend)
  end function integrate2D
end program
