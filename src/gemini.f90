Program Gemini3D
!! MAIN PROGRAM FOR GEMINI3D
!! Need program statement for FORD
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use sanity_check, only : check_finite_output
use phys_consts, only : lnchem, lwave, lsp, wp, debug
use grid, only: grid_size,read_grid,clear_grid, lx1,lx2,lx3,lx2all,lx3all
use mesh, only: curvmesh
use config, only : read_configfile, gemini_cfg, get_compiler_vendor
use pathlib, only : assert_file_exists, assert_directory_exists, expanduser, get_suffix
use io, only : input_plasma,create_outdir,output_plasma,create_outdir_aur,output_aur,find_milestone
use mpimod, only : mpisetup, mpibreakdown, mpi_manualgrid, mpigrid, lid, lid2,lid3,myid,myid2,myid3
use multifluid, only : fluid_adv
use neutral, only : neutral_atmos,make_dneu,neutral_perturb,clear_dneu,init_neutrals
use potentialBCs_mumps, only: clear_potential_fileinput, init_Efieldinput
use potential_comm,only : electrodynamics,pot2perpfield,velocities, get_BGEfields
use collisions, only: conductivities
use precipBCs_mod, only: clear_precip_fileinput, init_precipinput
use temporal, only : dt_comm
use timeutils, only: dateinc, find_lastdate

implicit none (type, external)

integer :: ierr

!> VARIABLES READ IN FROM CONFIG FILE
real(wp) :: UTsec
!! UT (s)
integer, dimension(3) :: ymd
!! year, month, day (current, not to be confused with starting year month and day in gemini_cfg structure)

type(gemini_cfg) :: cfg
!! holds many user simulation parameters

!> GRID STRUCTURE
type(curvmesh) :: x
!! structure containg grid locations, finite differences, etc.:  see grid module for details

!> STATE VARIABLES
!> MZ note:  it is likely that there could be a plasma and neutral derived type containing these data...  May be worth considering in a refactor...
real(wp), dimension(:,:,:,:), allocatable :: ns,vs1,vs2,vs3,Ts
!! fluid state variables
real(wp), dimension(:,:,:), allocatable :: E1,E2,E3,J1,J2,J3,Phi
!! electrodynamic state variables
real(wp), dimension(:,:,:), allocatable :: rhov2,rhov3,B1,B2,B3
!! inductive state vars. (for future use - except for B1 which is used for the background field)
real(wp), dimension(:,:,:), allocatable :: rhom,v1,v2,v3
!! inductive auxiliary
real(wp), dimension(:,:,:,:), allocatable :: nn
!! neutral density array
real(wp), dimension(:,:,:), allocatable :: Tn,vn1,vn2,vn3
!! neutral temperature and velocities
real(wp), dimension(:,:,:), allocatable :: Phiall
!! full-grid potential solution.  To store previous time step value
real(wp), dimension(:,:,:), allocatable :: iver
!! integrated volume emission rate of aurora calculated by GLOW

!TEMPORAL VARIABLES
real(wp) :: t=0, dt=1e-6_wp,dtprev
!! time from beginning of simulation (s) and time step (s)
real(wp) :: tout
!! time for next output and time between outputs
real(wp) :: tstart,tfin
!! temp. vars. for measuring performance of code blocks
integer :: it,isp, iupdate
!! time and species loop indices
real(wp) :: tneuBG !for testing whether we should re-evaluate neutral background

!> WORK ARRAYS
real(wp), allocatable :: dl1,dl2,dl3     !these are grid distances in [m] used to compute Courant numbers

real(wp) :: tglowout
!! time for next GLOW output

!> FOR HANDLING OUTPUT
integer :: lid2in,lid3in

!> TO CONTROL THROTTLING OF TIME STEP
real(wp), parameter :: dtscale=2.0_wp

!> Temporary variable for toggling full vs. other output
integer :: flagoutput
real(wp) :: tmilestone=0._wp

!> Milestone information
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp,ttmp,tdur
character(:), allocatable :: filetmp

!> For reproducing initial drifts; these are allocated and the deallocated since they can be large
real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
real(wp), dimension(:,:,:,:), allocatable :: muP,muH,muPvn,muHvn
real(wp), dimension(:,:,:), allocatable :: E01,E02,E03


!! MAIN PROGRAM
early_catch : block
integer i
character(512) :: argv
call get_command_argument(1, argv, status=i)
if (i<0) error stop 'command line truncated, plase file GitHub Issue'
if (i/=0 .or. argv == '-h' .or. argv == '-help') then
  print '(A,i4)', argv, i
  call help_cli()
endif
if (argv == '-compiler') then
  print '(A)', get_compiler_vendor()
  stop
endif
end block early_catch

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()

call cli(cfg, lid2in, lid3in)
!! initial_config is AFTER mpi_setup

!> CHECK THE GRID SIZE AND ESTABLISH A PROCESS GRID
call grid_size(cfg%indatsize)

!> MPI gridding cannot be done until we know the grid size
if (lid2in==-1) then
  call mpigrid(lx2all, lx3all)
  !! grid_size defines lx2all and lx3all
  print '(A, 2I6)', 'process grid (Number MPI processes) x2, x3:  ',lid2,lid3
  print '(A, I6, A, 2I6)', 'Process:',myid,' at process grid location:',myid2,myid3
else
  call mpi_manualgrid(lx2all, lx3all, lid2in, lid3in)
endif

!> LOAD UP THE GRID STRUCTURE/MODULE VARS. FOR THIS SIMULATION
call read_grid(cfg%indatsize,cfg%indatgrid,cfg%flagperiodic, x)
!! read in a previously generated grid from filenames listed in input file

!> CREATE/PREP OUTPUT DIRECTORY AND OUTPUT SIMULATION SIZE AND GRID DATA
!> ONLY THE ROOT PROCESS WRITES OUTPUT DATA

if (myid==0) then
  call create_outdir(cfg)
  if (cfg%flagglow /= 0) call create_outdir_aur(cfg%outdir)
end if


!> ALLOCATE ARRAYS (AT THIS POINT ALL SIZES ARE SET FOR EACH PROCESS SUBGRID)
allocate(ns(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs1(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs2(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), &
  vs3(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), Ts(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp))
allocate(rhov2(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhov3(-1:lx1+2,-1:lx2+2,-1:lx3+2),B1(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
         B2(-1:lx1+2,-1:lx2+2,-1:lx3+2),B3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(v1(-1:lx1+2,-1:lx2+2,-1:lx3+2),v2(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
         v3(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhom(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(E1(lx1,lx2,lx3),E2(lx1,lx2,lx3),E3(lx1,lx2,lx3),J1(lx1,lx2,lx3),J2(lx1,lx2,lx3),J3(lx1,lx2,lx3))
allocate(Phi(lx1,lx2,lx3))
allocate(nn(lx1,lx2,lx3,lnchem),Tn(lx1,lx2,lx3),vn1(lx1,lx2,lx3), vn2(lx1,lx2,lx3),vn3(lx1,lx2,lx3))


!> ALLOCATE MEMORY FOR ROOT TO STORE CERTAIN VARS. OVER ENTIRE GRID
if (myid==0) then
  allocate(Phiall(lx1,lx2all,lx3all))
end if

!> ALLOCATE MEMORY FOR AURORAL EMISSIONS, IF CALCULATED
if (cfg%flagglow /= 0) then
  allocate(iver(lx2,lx3,lwave))
  iver = 0
end if


!> Set initial time variables to simulation; this requires detecting whether we are trying to restart a simulation run
!> LOAD ICS AND DISTRIBUTE TO WORKERS (REQUIRES GRAVITY FOR INITIAL GUESSING)
!> ZZZ - this also should involve setting of Phiall...  Either to zero or what the input file specifies...
!        does not technically need to be broadcast to workers (since root sets up electrodynamics), but perhaps
!        should be anyway since that is what the user probably would expect and there is little performance penalty.
call find_milestone(cfg%outdir,get_suffix(cfg%indatsize),cfg%ymd0,cfg%UTsec0,cfg%dtout,ttmp,ymdtmp,UTsectmp,filetmp)
if (myid==0) print*, 'Last milestone (if any)found in output directory:  ',ymdtmp,UTsectmp,filetmp
if ( any(ymdtmp/=cfg%ymd0) .or. abs(UTsectmp-cfg%UTsec0)>cfg%dtout ) then  !! treat this as a restart scenario
  if (myid==0) then
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*, '! Restarting simulation from time:  ',ymdtmp,UTsectmp
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if

  !! Set start variables accordingly and read in the milestone
  UTsec=UTsectmp
  ymd=ymdtmp
  tdur=cfg%tdur-ttmp    ! subtract off time that has elapsed to milestone
  if (myid==0) then
    print*, 'Treating the following file as initial conditions:  ',filetmp
    print*, ' full duration:  ',cfg%tdur,'; remaining simulation time:  ',tdur
  end if
  cfg%tdur=tdur         ! just to insure consistency
  call input_plasma(x%x1,x%x2all,x%x3all,cfg%indatsize,filetmp,ns,vs1,Ts,Phi,Phiall)
else !! start at the beginning
  UTsec = cfg%UTsec0
  ymd = cfg%ymd0
  tdur = cfg%tdur
  call input_plasma(x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,Phiall)
end if
it = 1
t = 0
tout = t
tglowout = t
tneuBG=t


!ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
!THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
!FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INTERTIAL CAPACITANCE THERE
!WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
!BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
!VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
rhov2 = 0
rhov3 = 0
v2 = 0
v3 = 0
B2 = 0
B3 = 0
B1(1:lx1,1:lx2,1:lx3) = x%Bmag
!! this assumes that the grid is defined s.t. the x1 direction corresponds
!! to the magnetic field direction (hence zero B2 and B3).


!> Inialize neutral atmosphere, note the use of fortran's weird scoping rules to avoid input args.  Must occur after initial time info setup
if(myid==0) print*, 'Priming neutral input'
call init_neutrals(dt,t,cfg,ymd,UTsec,x,nn,Tn,vn1,vn2,vn3)

!> Initialize auroral inputs; must occur after initial timing info setup
if(myid==0) print*, 'Priming electric field input'
call init_Efieldinput(dt,t,cfg,ymd,UTsec,x)

if(myid==0) print*, 'Priming precipitation input'
call init_precipinput(dt,t,cfg,ymd,UTsec,x)


!> Recompute electrodynamic quantities needed for restarting
!! these do not include background
E1 = 0
call pot2perpfield(Phi,x,E2,E3)
if(myid==0) then
  print*, 'Recomputed initial dist. fields:  '
  print*, '    ',minval(E1),maxval(E1)
  print*, '    ',minval(E2),maxval(E2)
  print*, '    ',minval(E3),maxval(E3)
end if
allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))
E01=0; E02=0; E03=0;
if (cfg%flagE0file==1) then
  call get_BGEfields(x,E01,E02,E03)
end if
if(myid==0) then
  print*, 'Recomputed initial BG fields:  '
  print*, '    ',minval(E01),maxval(E01)
  print*, '    ',minval(E02),maxval(E02)
  print*, '    ',minval(E03),maxval(E03)
end if

allocate(sig0(lx1,lx2,lx3),sigP(lx1,lx2,lx3),sigH(lx1,lx2,lx3),sigPgrav(lx1,lx2,lx3),sigHgrav(lx1,lx2,lx3))
allocate(muP(lx1,lx2,lx3,lsp),muH(lx1,lx2,lx3,lsp),muPvn(lx1,lx2,lx3,lsp),muHvn(lx1,lx2,lx3,lsp))
call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH,muPvn,muHvn,sigPgrav,sigHgrav)
E1=E1+E01; E2=E2+E02; E3=E3+E03
call velocities(muP,muH,muPvn,muHvn,E2,E3,vn2,vn3,cfg%flaggravdrift,vs2,vs3)
deallocate(sig0,sigP,sigH,muP,muH,muPvn,muHvn,sigPgrav,sigHgrav)
deallocate(E01,E02,E03)
if(myid==0) then
  print*, 'Recomputed initial drifts:  '
  print*, '    ',minval(vs2),maxval(vs2)
  print*, '    ',minval(vs3),maxval(vs3)
end if

!> control update rate from excessive console printing
!! considering small vs. large simulations
!! these are arbitrary levels, so feel free to finesse
if (lx1*lx2*lx3 < 20000) then
  iupdate = 50
elseif (lx1*lx2*lx3 < 100000) then
  iupdate = 10
else
  iupdate = 1
endif


!> Main time loop
do while (t < tdur)
  !! TIME STEP CALCULATION, requires workers to report their most stringent local stability constraint
  dtprev = dt
  call dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)
  if (it>1) then
    if(dt/dtprev > dtscale) then
      !! throttle how quickly we allow dt to increase
      dt=dtscale*dtprev
      if (myid == 0) then
        print '(A,EN14.3)', 'Throttling dt to:  ',dt
      end if
    end if
  end if

  !COMPUTE BACKGROUND NEUTRAL ATMOSPHERE USING MSIS00.
  if ( it/=1 .and. cfg%flagneuBG .and. t>tneuBG) then     !we dont' throttle for tneuBG so we have to do things this way to not skip over...
    call cpu_time(tstart)
    call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,nn,Tn,vn1,vn2,vn3)
    tneuBG=tneuBG+cfg%dtneuBG;
    if (myid==0) then
      call cpu_time(tfin)
      print *, 'Neutral background at time:  ',t,' calculated in time:  ',tfin-tstart
    end if
  end if


  !> GET NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
  if (cfg%flagdneu==1) then
    call cpu_time(tstart)
    call neutral_perturb(cfg,dt,cfg%dtneu,t,ymd,UTsec,x,nn,Tn,vn1,vn2,vn3)
    if (myid==0 .and. debug) then
      call cpu_time(tfin)
      print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
    endif
  end if

  !! POTENTIAL SOLUTION
  call cpu_time(tstart)
  call electrodynamics(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x,E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
  if (myid==0 .and. debug) then
    call cpu_time(tfin)
    print *, 'Electrodynamics total solve time:  ',tfin-tstart
  endif

  !> UPDATE THE FLUID VARIABLES
  if (myid==0 .and. debug) call cpu_time(tstart)
  call fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,cfg,t,dt,x,nn,vn1,vn2,vn3,Tn,iver,ymd,UTsec)
  if (myid==0 .and. debug) then
    call cpu_time(tfin)
    print *, 'Multifluid total solve time:  ',tfin-tstart
  endif

  !> FIXME:  MZ - shouldn't this be done for all workers; also how much overhead does this incur every time step???
  !> Sanity check key variables before advancing
  if (myid==0) call check_finite_output(t, myid, vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)

  !> NOW OUR SOLUTION IS FULLY UPDATED SO UPDATE TIME VARIABLES TO MATCH...
  it = it + 1
  t = t + dt
  if (myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',cfg%tdur
  call dateinc(dt,ymd,UTsec)

  if (myid==0 .and. (modulo(it, iupdate) == 0 .or. debug)) then
    !! print every 10th time step to avoid extreme amounts of console printing
    print '(A,I4,A1,I0.2,A1,I0.2,A1,F12.6,A5,F8.6)', 'Current time ',ymd(1),'-',ymd(2),'-',ymd(3),' ',UTsec,'; dt=',dt
  endif

  if (cfg%dryrun) then
    ierr = mpibreakdown()
    if (ierr /= 0) error stop 'Gemini dry run MPI shutdown failure'
    block
      character(8) :: date
      character(10) :: time

      call date_and_time(date,time)
      print '(/,A)', 'DONE: ' // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' &
        // time(1:2) // ':' // time(3:4) // ':' // time(5:)
      stop "OK: Gemini dry run"
    end block
  endif

  !> File output
  if (abs(t-tout) < 1d-5) then
    tout = tout + cfg%dtout
    if (cfg%nooutput .and. myid==0) then
      write(stderr,*) 'WARNING: skipping file output at sim time (sec)',t
      cycle
    endif
    !! close enough to warrant an output now...
    if (myid==0 .and. debug) call cpu_time(tstart)

    !! We may need to adjust flagoutput if we are hitting a milestone
    flagoutput=cfg%flagoutput
    if (cfg%mcadence>0 .and. abs(t-tmilestone) < 1d-5) then
      flagoutput=1    !force a full output at the milestone
      call output_plasma(cfg%outdir,flagoutput,ymd, &
        UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
        cfg%out_format)
      tmilestone=t+cfg%dtout*cfg%mcadence
      if(myid==0) print*, 'Milestone output triggered.'
    else
      call output_plasma(cfg%outdir,flagoutput,ymd, &
        UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
        cfg%out_format)
    end if
    if (myid==0 .and. debug) then
      call cpu_time(tfin)
      print *, 'Plasma output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    endif
  end if

  !> GLOW file output
  if ((cfg%flagglow /= 0) .and. (abs(t-tglowout) < 1d-5)) then !same as plasma output
    call cpu_time(tstart)
    call output_aur(cfg%outdir, cfg%flagglow, ymd, UTsec, iver, cfg%out_format)
    if (myid==0) then
      call cpu_time(tfin)
      print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
    end if
    tglowout = tglowout + cfg%dtglowout
  end if
end do


!! DEALLOCATE MAIN PROGRAM DATA
deallocate(ns,vs1,vs2,vs3,Ts)
deallocate(E1,E2,E3,J1,J2,J3)
deallocate(nn,Tn,vn1,vn2,vn3)

if (myid==0) deallocate(Phiall)

if (cfg%flagglow/=0) deallocate(iver)

!! DEALLOCATE MODULE VARIABLES (MAY HAPPEN AUTOMATICALLY IN F2003???)
call clear_grid(x)
call clear_dneu()
call clear_precip_fileinput()
call clear_potential_fileinput()
!call clear_BGfield()


!! SHUT DOWN MPI
ierr = mpibreakdown()

if (ierr /= 0) then
  write(stderr, *) 'GEMINI: abnormal MPI shutdown code', ierr, 'Process #', myid,' /',lid-1
  error stop
endif

block
  character(8) :: date
  character(10) :: time

  call date_and_time(date,time)
  print '(/,A,I6,A,I6,A)', 'GEMINI normal termination, Process #', myid,' /',lid-1, ' at ' // date // 'T' // time
end block

contains


subroutine cli(cfg, lid2in, lid3in)

type(gemini_cfg), intent(out) :: cfg
integer, intent(out) :: lid2in, lid3in

integer :: argc, i
character(512) :: argv
character(8) :: date
character(10) :: time

argc = command_argument_count()


if(lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'

call get_command_argument(0, argv)
call date_and_time(date,time)
print '(2A,I6,A3,I6,A)', trim(argv), ' Process:  ', myid,' / ',lid-1, ' at ' // date // 'T' // time


!> READ FILE INPUT
call get_command_argument(1, argv, status=i)
if (i/=0) error stop 'bad command line'
cfg%outdir = expanduser(trim(argv))

find_cfg : block
logical :: exists
character(*), parameter :: locs(4) = [character(18) :: "/inputs/config.nml", "/config.nml", "/inputs/config.ini", "/config.ini"]
character(:), allocatable :: loc
do i = 1,size(locs)
  loc = trim(cfg%outdir // locs(i))
  inquire(file=loc, exist=exists)
  if (exists) then
    cfg%infile = loc
    exit find_cfg
  endif
end do

write(stderr,*) 'gemini.bin: could not find config file in ',cfg%outdir
error stop 6
end block find_cfg

call read_configfile(cfg, verbose=.false.)

!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (myid==0) then
  call assert_file_exists(cfg%indatsize)
  call assert_file_exists(cfg%indatgrid)
  call assert_file_exists(cfg%indatfile)

  print *, '******************** input config ****************'
  print '(A)', 'simulation directory: ' // cfg%outdir
  print '(A51,I6,A1,I0.2,A1,I0.2)', ' start year-month-day:  ', cfg%ymd0(1), '-', cfg%ymd0(2),'-', cfg%ymd0(3)
  print '(A51,F10.3)', 'start time:  ',cfg%UTsec0
  print '(A51,F10.3)', 'duration:  ',cfg%tdur
  print '(A51,F10.3)', 'output every:  ',cfg%dtout
  print '(A,/,A,/,A,/,A)', 'gemini.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

  if(cfg%flagdneu==1) then
    call assert_directory_exists(cfg%sourcedir)
    print *, 'Neutral disturbance mlat,mlon:  ',cfg%sourcemlat,cfg%sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',cfg%dtneu
    print *, 'Neutral grid resolution (m):  ',cfg%drhon,cfg%dzn
    print *, 'Neutral disturbance data files located in directory:  ',cfg%sourcedir
  else
    print *, "no neutral disturbance specified."
  end if

  if (cfg%flagprecfile==1) then
    call assert_directory_exists(cfg%precdir)
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',cfg%dtprec
    print *, 'Precipitation file input source directory:  ' // cfg%precdir
  else
    print *, "no precipitation specified"
  end if

  if(cfg%flagE0file==1) then
    call assert_directory_exists(cfg%E0dir)
    print *, 'Electric field file input cadence (s):  ',cfg%dtE0
    print *, 'Electric field file input source directory:  ' // cfg%E0dir
  else
    print *, "no Efield specified"
  end if

  if (cfg%flagglow==1) then
    print *, 'GLOW enabled for auroral emission calculations.'
    print *, 'GLOW electron transport calculation cadence (s): ', cfg%dtglow
    print *, 'GLOW auroral emission output cadence (s): ', cfg%dtglowout
  else
    print *, "GLOW disabled"
  end if

  if (cfg%flagEIA) then
    print*, 'EIA enables with peok equatorial drift:  ',cfg%v0equator
  else
    print*, 'EIA disabled'
  end if

  if (cfg%flagneuBG) then
    print*, 'Variable background neutral atmosphere enabled at cadence:  ',cfg%dtneuBG
  else
    print*, 'Variable background neutral atmosphere disabled.'
  end if

  print*, 'Background precipitation has total energy flux and energy:  ',cfg%PhiWBG,cfg%W0BG

  if (cfg%flagJpar) then
    print*, 'Parallel current calculation enabled.'
  else
    print*, 'Parallel current calculation disabled.'
  end if

  print*, 'Inertial capacitance calculation type:  ',cfg%flagcap

  print*, 'Diffusion solve type:  ',cfg%diffsolvetype

  if (cfg%mcadence > 0) then
    print*, 'Milestone output selected; cadence (every nth outout) of:  ',cfg%mcadence
  else
    print*, 'Milestone output disabled.'
  end if

  if (cfg%flaggravdrift) then
    print*, 'Gravitational drift terms enabled.'
  else
    print*, 'Gravitaional drift terms disabled.'
  end if

  print *,  '**************** end input config ***************'
end if

!! default values
lid2in = -1  !< sentinel

do i = 2,argc
  call get_command_argument(i,argv)

  select case (argv)
  case ('-h', '-help')
    call help_cli()
  case ('-compiler')
    print '(A)', get_compiler_vendor()
    stop
  case ('-d', '-debug')
    debug = .true.
  case ('-dryrun')
    !! this is a no file output test mode that runs one time step then quits
    !! it helps avoid HPC queuing when a simple setup error exists
    cfg%dryrun = .true.
  case ('-nooutput')
    cfg%nooutput = .true.
  case ('-out_format')
    !! used mostly for debugging--normally should be set as file_format in config.nml
    call get_command_argument(i+1, argv, status=ierr)
    if(ierr/=0) error stop 'gemini.bin -out_format {h5,nc,dat} parameter is required'
    cfg%out_format = trim(argv)
    print *,'override output file format: ',cfg%out_format
  case ('-manual_grid')
    call get_command_argument(i+1, argv, status=ierr)
    if(ierr/=0) error stop 'gemini.bin -manual_grid lx2 lx3 parameters are required'
    read(argv,*) lid2in
    call get_command_argument(i+2, argv, status=ierr)
    if(ierr/=0) error stop 'gemini.bin -manual_grid lx2 lx3 parameters are required'
    read(argv,*) lid3in
  end select
end do

end subroutine cli

subroutine help_cli()
print '(/,A,/)', 'GEMINI-3D: by Matthew Zettergren'
print '(A)', 'GLOW and auroral interfaces by Guy Grubbs'
print '(A)', 'build system and software engineering by Michael Hirsch'
print '(A,/)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'must specify simulation output directory. Example:'
print '(/,A,/)', '  mpiexec -np 4 build/gemini.bin /path/to/simulation_outputs'
print '(A)', '-dryrun    allows quick check of first time step'
print '(A)', '-manual_grid lx2 lx3    defines the number of MPI processes along x2 and x3.'
print '(A)', '  If -manual_grid is not specified, the MPI processes are auto-assigned along x2 and x3.'
stop 'EOF: Gemini-3D'
end subroutine help_cli

end program
