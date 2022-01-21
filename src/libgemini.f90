! Copyright 2021 Matthew Zettergren

! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!   http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module gemini3d
! top-level module for Gemini3D
use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_bool, c_float
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use gemini_init, only : find_config, check_input_files
use gemini_cli, only : cli
use config, only : read_configfile
use sanity_check, only : check_finite_output, check_finite_pertub
use phys_consts, only : lnchem, lwave, lsp, wp, debug
use grid, only: grid_size,read_grid,grid_drift, lx1,lx2,lx3,lx2all,lx3all
use meshobj, only: curvmesh
use config, only : gemini_cfg
use io, only : input_plasma,create_outdir,output_plasma, output_aur, output_cond, find_milestone
use filesystem, only : expanduser
use mpimod, only : mpisetup, mpibreakdown, mpi_manualgrid, process_grid_auto, mpi_cfg
use multifluid, only : fluid_adv

use msis_interface, only : msisinit
use neutral, only : neutral_atmos,make_dneu,neutral_perturb,clear_dneu,init_neutrals, neutral_winds

use potentialBCs_mumps, only: init_Efieldinput
use potential_comm,only : electrodynamics,pot2perpfield,velocities, get_BGEfields
use collisions, only: conductivities
use precipBCs_mod, only: init_precipinput
use temporal, only : dt_comm
use timeutils, only: dateinc, find_lastdate

implicit none (type, external)

type, bind(C) :: c_params
!! this MUST match gemini3d.h and gemini_main.f90 exactly including order
logical(c_bool) :: fortran_cli, debug, dryrun
character(kind=c_char) :: out_dir(1000)
!! .ini [base]
integer(c_int) :: ymd(3)
real(kind=c_float) :: UTsec0, tdur, dtout, activ(3), tcfl, Teinf
!! .ini
end type c_params

contains


subroutine gemini_main(p, lid2in, lid3in)  bind(C)
!! NOTE: if use_cli=.true., then {out_dir, lid2in, lid3in} are ignored and CLI is used instead.

type(c_params), intent(in) :: p
!! output directory for Gemini3D to write simulation data to (can be large files GB, TB, ...)
integer(c_int), intent(inout) :: lid2in, lid3in  !< inout to allow optional CLI

integer :: ierr
logical :: exists

!> VARIABLES READ IN FROM CONFIG FILE
real(wp) :: UTsec
!! UT (s)
integer, dimension(3) :: ymd
!! year, month, day (current, not to be confused with starting year month and day in gemini_cfg structure)

type(gemini_cfg) :: cfg
!! holds many user simulation parameters

!> grid type (polymorphic) containing geometric information and associate procedures
class(curvmesh), allocatable :: x

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

!> TO CONTROL THROTTLING OF TIME STEP
real(wp), parameter :: dtscale=2

!> Temporary variable for toggling full vs. other output
integer :: flagoutput
real(wp) :: tmilestone = 0

!> Milestone information
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp,ttmp,tdur
character(:), allocatable :: filetmp

!> For reproducing initial drifts; these are allocated and the deallocated since they can be large
real(wp), dimension(:,:,:), allocatable :: sig0,sigP,sigH,sigPgrav,sigHgrav
real(wp), dimension(:,:,:,:), allocatable :: muP,muH,nusn
real(wp), dimension(:,:,:), allocatable :: E01,E02,E03

!> Describing Lagrangian grid (if used)
real(wp) :: v2grid,v3grid

character(*), parameter :: msis2_param_file = "msis20.parm"

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()

if(mpi_cfg%lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'

if(p%fortran_cli) then
  call cli(cfg, lid2in, lid3in, debug)
else
  block
    character(size(p%out_dir)) :: buf
    integer :: i
    buf = "" !< ensure buf has no garbage characters

    do i = 1, len(buf)
      if (p%out_dir(i) == c_null_char) exit
      buf(i:i) = p%out_dir(i)
    enddo
    cfg%outdir = expanduser(buf)

    cfg%dryrun = p%dryrun
    debug = p%debug
  end block
endif

call find_config(cfg)

call read_configfile(cfg, verbose=.false.)

call check_input_files(cfg)


!> CHECK THE GRID SIZE AND ESTABLISH A PROCESS GRID
call grid_size(cfg%indatsize)

!> MPI gridding cannot be done until we know the grid size
if (lid2in==-1) then
  call process_grid_auto(lx2all, lx3all)
  !! grid_size defines lx2all and lx3all
else
  call mpi_manualgrid(lx2all, lx3all, lid2in, lid3in)
endif
print '(A, I0, A1, I0)', 'process grid (Number MPI processes) x2, x3:  ',mpi_cfg%lid2, ' ', mpi_cfg%lid3
print '(A, I0, A, I0, A1, I0)', 'Process:',mpi_cfg%myid,' at process grid location: ',mpi_cfg%myid2,' ',mpi_cfg%myid3

!> LOAD UP THE GRID STRUCTURE/MODULE VARS. FOR THIS SIMULATION
call read_grid(cfg%indatsize,cfg%indatgrid,cfg%flagperiodic, x)
!! read in a previously generated grid from filenames listed in input file

!> CREATE/PREP OUTPUT DIRECTORY AND OUTPUT SIMULATION SIZE AND GRID DATA
!> ONLY THE ROOT PROCESS WRITES OUTPUT DATA

if (mpi_cfg%myid==0) then
  call create_outdir(cfg)
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
if (mpi_cfg%myid==0) then
  allocate(Phiall(lx1,lx2all,lx3all))
end if

!> ALLOCATE MEMORY FOR AURORAL EMISSIONS, IF CALCULATED
if (cfg%flagglow /= 0) then
  allocate(iver(lx2,lx3,lwave))
  iver = 0
end if


!> FIXME: Zero out all state variables here - the corner ghost cells otherwise never get set and could contain garbage whicgh may make the sanity check fail since it does look at ghost cells, as well.
!ns=0; vs1=0; Ts=0;
! oddly this doesn't seem to help our issue...

!> Set initial time variables to simulation; this requires detecting whether we are trying to restart a simulation run
!> LOAD ICS AND DISTRIBUTE TO WORKERS (REQUIRES GRAVITY FOR INITIAL GUESSING)
!> ZZZ - this also should involve setting of Phiall...  Either to zero or what the input file specifies...
!        does not technically need to be broadcast to workers (since root sets up electrodynamics), but perhaps
!        should be anyway since that is what the user probably would expect and there is little performance penalty.
call find_milestone(cfg, ttmp, ymdtmp, UTsectmp, filetmp)
if ( ttmp > 0 ) then
  !! restart scenario
  if (mpi_cfg%myid==0) then
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print*, '! Restarting simulation from time:  ',ymdtmp,UTsectmp
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if

  !! Set start variables accordingly and read in the milestone
  UTsec=UTsectmp
  ymd=ymdtmp
  tdur=cfg%tdur-ttmp    ! subtract off time that has elapsed to milestone
  if (mpi_cfg%myid==0) then
    print*, 'Treating the following file as initial conditions:  ',filetmp
    print*, ' full duration:  ',cfg%tdur,'; remaining simulation time:  ',tdur
  end if

  if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Cannot restart simulation from the final time step!'

  cfg%tdur=tdur         ! just to insure consistency

  call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,filetmp,ns,vs1,Ts,Phi,Phiall)
else !! start at the beginning
  UTsec = cfg%UTsec0
  ymd = cfg%ymd0
  tdur = cfg%tdur

  if (tdur <= 1e-6_wp .and. mpi_cfg%myid==0) error stop 'Simulation is of zero time duration'

  call input_plasma(cfg%outdir, x%x1,x%x2all,x%x3all,cfg%indatsize,cfg%indatfile,ns,vs1,Ts,Phi,Phiall)
end if

it = 1
t = 0
tout = t
tglowout = t
tneuBG=t


!ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
!THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
!FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INERTIAL CAPACITANCE THERE
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
if(mpi_cfg%myid==0) print*, 'Priming electric field input'
call init_Efieldinput(dt,t,cfg,ymd,UTsec,x)

allocate(E01(lx1,lx2,lx3),E02(lx1,lx2,lx3),E03(lx1,lx2,lx3))
E01=0; E02=0; E03=0;
if (cfg%flagE0file==1) then
  call get_BGEfields(x,E01,E02,E03)
end if
if (cfg%flaglagrangian) then    ! Lagrangian (moving) grid; compute from input background electric fields
  call grid_drift(x,E02,E03,v2grid,v3grid)
  if (mpi_cfg%myid==0) print*, mpi_cfg%myid,' using Lagrangian grid moving at:  ',v2grid,v3grid
else                            ! stationary grid
  v2grid = 0
  v3grid = 0
  E1 = E1 + E01
  E2 = E2 + E02
  E3 = E3 + E03
end if

!> Precipitation input setup
if(mpi_cfg%myid==0) print*, 'Priming precipitation input'
call init_precipinput(dt,t,cfg,ymd,UTsec,x)

!> Neutral atmosphere setup
if(cfg%msis_version == 20) then
  inquire(file=msis2_param_file, exist=exists)
  if(.not.exists) error stop 'could not find MSIS 2 parameter file ' // msis2_param_file // &
    ' this file must be in the same directory as gemini.bin, and run from that directory. ' // &
    'This limitation comes from how MSIS 2.x is coded internally.'
  call msisinit(parmfile=msis2_param_file)
end if

if(mpi_cfg%myid==0) print*, 'Computing background and priming neutral perturbation input (if used)'
call init_neutrals(dt,t,cfg,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)

!> Recompute electrodynamic quantities needed for restarting
!> these do not include background
E1 = 0
call pot2perpfield(Phi,x,E2,E3)
if(mpi_cfg%myid==0) then
  print '(A)', 'Recomputed initial dist. fields:'
  print*, '    gemini ',minval(E1),maxval(E1)
  print*, '    gemini ',minval(E2),maxval(E2)
  print*, '    gemini ',minval(E3),maxval(E3)

  print*, 'Recomputed initial BG fields:'
  print*, '    ',minval(E01),maxval(E01)
  print*, '    ',minval(E02),maxval(E02)
  print*, '    ',minval(E03),maxval(E03)
end if


!> Recompute drifts and make some decisions about whether to invoke a Lagrangian grid
allocate(sig0(lx1,lx2,lx3))
allocate(sigP, sigH, sigPgrav, sigHgrav, mold=sig0)
allocate(muP(lx1,lx2,lx3,lsp))
allocate(muH, nusn, mold=muP)
call conductivities(nn,Tn,ns,Ts,vs1,B1, sig0,sigP,sigH,muP,muH,nusn,sigPgrav,sigHgrav)
call velocities(muP,muH,nusn,E2,E3,vn2,vn3,ns,Ts,x,cfg%flaggravdrift,cfg%flagdiamagnetic,vs2,vs3)
!> optional output to file
call output_cond(cfg%outdir, ymd, UTsec, sig0, sigP, sigH, cfg%out_format)
deallocate(sig0, sigP, sigH, muP, muH, nusn, sigPgrav, sigHgrav)
deallocate(E01,E02,E03)
if(mpi_cfg%myid == 0) then
  print*, 'Recomputed initial drifts:  '
  print*, '    ',minval(vs2(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs2(1:lx1,1:lx2,1:lx3,1:lsp))
  print*, '    ',minval(vs3(1:lx1,1:lx2,1:lx3,1:lsp)),maxval(vs3(1:lx1,1:lx2,1:lx3,1:lsp))
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
main : do while (t < tdur)
  !> TIME STEP CALCULATION, requires workers to report their most stringent local stability constraint
  dtprev = dt
  call dt_comm(t,tout,tglowout,cfg,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,dt)
  if (it>1) then
    if(dt/dtprev > dtscale) then
      !! throttle how quickly we allow dt to increase
      dt=dtscale*dtprev
      if (mpi_cfg%myid == 0) then
        print '(A,EN14.3)', 'Throttling dt to:  ',dt
      end if
    end if
  end if

  !> COMPUTE BACKGROUND NEUTRAL ATMOSPHERE USING MSIS
  if ( it/=1 .and. cfg%flagneuBG .and. t>tneuBG) then     !we dont' throttle for tneuBG so we have to do things this way to not skip over...
    call cpu_time(tstart)
    call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,cfg%activ,nn,Tn,cfg%msis_version)
    call neutral_winds(ymd, UTsec, Ap=cfg%activ(3), x=x, v2grid=v2grid,v3grid=v3grid,vn1=vn1, vn2=vn2, vn3=vn3)
    tneuBG=tneuBG+cfg%dtneuBG;
    if (mpi_cfg%myid==0) then
      call cpu_time(tfin)
      print *, 'Neutral background at time:  ',t,' calculated in time:  ',tfin-tstart
    end if
  end if

  !> GET NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
  if (cfg%flagdneu==1) then
    call cpu_time(tstart)
    call neutral_perturb(cfg,dt,cfg%dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)
    call check_finite_pertub(cfg%outdir, t, mpi_cfg%myid, nn, Tn, vn1, vn2, vn3)

    if (mpi_cfg%myid==0 .and. debug) then
      call cpu_time(tfin)
      print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
    endif
  end if


  !> POTENTIAL SOLUTION
  call cpu_time(tstart)
  call electrodynamics(it,t,dt,nn,vn2,vn3,Tn,cfg,ns,Ts,vs1,B1,vs2,vs3,x,E1,E2,E3,J1,J2,J3,Phiall,ymd,UTsec)
  if (mpi_cfg%myid==0 .and. debug) then
    call cpu_time(tfin)
    print *, 'Electrodynamics total solve time:  ',tfin-tstart
  endif

  !> UPDATE THE FLUID VARIABLES
  if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)
  call fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,cfg,t,dt,x,nn,vn1,vn2,vn3,Tn,iver,ymd,UTsec, first=(it==1) )
  if (mpi_cfg%myid==0 .and. debug) then
    call cpu_time(tfin)
    print *, 'Multifluid total solve time:  ',tfin-tstart
  endif

  !> Sanity check key variables before advancing
  ! FIXME: for whatever reason, it is just a fact that vs1 has trash in ghost cells after fluid_adv; I don't know why...
  call check_finite_output(cfg%outdir, t, mpi_cfg%myid, vs2,vs3,ns,vs1,Ts, Phi,J1,J2,J3)

  !> NOW OUR SOLUTION IS FULLY UPDATED SO UPDATE TIME VARIABLES TO MATCH...
  it = it + 1
  t = t + dt
  if (mpi_cfg%myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',cfg%tdur
  call dateinc(dt,ymd,UTsec)

  if (mpi_cfg%myid==0 .and. (modulo(it, iupdate) == 0 .or. debug)) then
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
    if (cfg%nooutput ) then
      if (mpi_cfg%myid==0) write(stderr,*) 'WARNING: skipping file output at sim time (sec)',t
      cycle main
    endif
    !! close enough to warrant an output now...
    if (mpi_cfg%myid==0 .and. debug) call cpu_time(tstart)

    !! We may need to adjust flagoutput if we are hitting a milestone
    flagoutput=cfg%flagoutput
    if (cfg%mcadence>0 .and. abs(t-tmilestone) < 1d-5) then
      flagoutput=1    !force a full output at the milestone
      call output_plasma(cfg%outdir,flagoutput,ymd, &
        UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
        cfg%out_format)
      tmilestone = t + cfg%dtout * cfg%mcadence
      if(mpi_cfg%myid==0) print*, 'Milestone output triggered.'
    else
      call output_plasma(cfg%outdir,flagoutput,ymd, &
        UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, &
        cfg%out_format)
    end if
    if (mpi_cfg%myid==0 .and. debug) then
      call cpu_time(tfin)
      print *, 'Plasma output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    endif
  end if

  !> GLOW file output
  if ((cfg%flagglow /= 0) .and. (abs(t-tglowout) < 1d-5)) then !same as plasma output
    call cpu_time(tstart)
    call output_aur(cfg%outdir, cfg%flagglow, ymd, UTsec, iver, cfg%out_format)
    if (mpi_cfg%myid==0) then
      call cpu_time(tfin)
      print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
    end if
    tglowout = tglowout + cfg%dtglowout
  end if
end do main


!> DEALLOCATE module data. We haven't verified it's strictly necessary, but it's been our practice.
deallocate(ns,vs1,vs2,vs3,Ts)
deallocate(E1,E2,E3,J1,J2,J3)
deallocate(nn,Tn,vn1,vn2,vn3)

if (mpi_cfg%myid==0) deallocate(Phiall)

if (cfg%flagglow/=0) deallocate(iver)

!> DEALLOCATE MODULE VARIABLES (MAY HAPPEN AUTOMATICALLY IN F2003???)
!call clear_grid(x)
call clear_dneu()
! now taken care of by an object destructor...
!call clear_precip_fileinput()
!call clear_potential_fileinput()
!call clear_BGfield()

end subroutine gemini_main


end module gemini3d
