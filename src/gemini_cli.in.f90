module gemini_cli

use config, only : read_configfile, gemini_cfg, get_compiler_vendor
use pathlib, only : assert_file_exists, assert_directory_exists, expanduser
use mpimod, only : mpisetup, mpibreakdown

implicit none (type, external)
private
public :: cli

contains


subroutine cli(myid, lid, cfg, lid2in, lid3in, debug)

integer, intent(in) :: myid, lid
type(gemini_cfg), intent(out) :: cfg
integer, intent(out) :: lid2in, lid3in
logical, intent(inout) :: debug

integer :: argc, i, ierr
character(512) :: argv
character(8) :: date
character(10) :: time

cfg%git_revision = "@git_rev@"

argc = command_argument_count()

call get_command_argument(1, argv, status=i)
if (i/=0) call help_cli(cfg%git_revision)

select case (argv)
case ('-h', '-help')
  call help_cli(cfg%git_revision)
case ('-compiler')
  print '(A)', get_compiler_vendor()
  stop
case ('-git')
  print '(A)', cfg%git_revision
  stop
end select

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()

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

error stop 'gemini.bin: could not find config file in ' // cfg%outdir

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

  if (cfg%flaglagrangian) then
    print*, 'Lagrangian grid enabled.'
  else
    print*, 'Lagrangian grid disabled'
  end if

  print *,  '**************** end input config ***************'
end if

!! default values
lid2in = -1  !< sentinel

do i = 2,argc
  call get_command_argument(i,argv)

  select case (argv)
  case ('-h', '-help')
    ierr = mpibreakdown()
    if (myid == 0) call help_cli(cfg%git_revision)
    stop
  case ('-compiler')
    ierr = mpibreakdown()
    if (myid==0) print '(A)', get_compiler_vendor()
    stop
  case ('-git')
    ierr = mpibreakdown()
    if (myid==0) print '(A)', cfg%git_revision
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

subroutine help_cli(git_revision)
character(*), intent(in) :: git_revision

print '(/,A,/)', 'GEMINI-3D ' // git_revision
print '(A)', 'by Matthew Zettergren'
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

end module gemini_cli
