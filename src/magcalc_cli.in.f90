module magcalc_cli

use, intrinsic :: iso_fortran_env, only : compiler_version
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
character(256) :: argv
character(8) :: date
character(10) :: time

logical :: file_exists

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


!> READ CONFIG FILE FROM OUTPUT DIRECTORY
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

if (cfg%outdir(1:1) == "-") then
  error stop 'magcalc.bin: not a known CLI option: ' // cfg%outdir
else
  error stop 'magcalc.bin: could not find config file in ' // cfg%outdir
endif

end block find_cfg
!> GRAB THE INFO FOR WHERE THE OUTPUT CALCULATIONS ARE STORED
call get_command_argument(2,argv)
cfg%fieldpointfile = trim(argv)
!! this file contains the field points at which we are computing magnetic perturbations, it will be copied into the output directory

if (myid==0) then
  print *, 'Simulation data directory:  ', cfg%outdir
  print *, 'Input config file:  ',cfg%infile
  print *, 'fieldpoint file: ', cfg%fieldpointfile
end if

call read_configfile(cfg)


!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (myid==0) then
  call assert_file_exists(cfg%indatsize)
  call assert_file_exists(cfg%indatgrid)
  call assert_file_exists(cfg%indatfile)

  print '(A,I6,A1,I0.2,A1,I0.2)', cfg%infile // ' start year-month-day:  ',cfg%ymd0(1),'-',cfg%ymd0(2),'-',cfg%ymd0(3)
  print '(A51,F10.3)', 'start time:  ',cfg%UTsec0
  print '(A51,F10.3)', 'duration:  ',cfg%tdur
  print '(A51,F10.3)', 'output every:  ',cfg%dtout
  print '(A,/,A,/,A,/,A)', 'magcalc.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

  if(cfg%flagdneu==1) then
    call assert_directory_exists(cfg%sourcedir)
    print *, 'Neutral disturbance mlat,mlon:  ',cfg%sourcemlat,cfg%sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',cfg%dtneu
    print *, 'Neutral grid resolution (m):  ',cfg%drhon,cfg%dzn
    print *, 'Neutral disturbance data files located in directory:  ',cfg%sourcedir
  end if

  if (cfg%flagprecfile==1) then
    call assert_directory_exists(cfg%precdir)
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',cfg%dtprec
    print *, 'Precipitation file input source directory:  ' // cfg%precdir
  end if

  if(cfg%flagE0file==1) then
    call assert_directory_exists(cfg%E0dir)
    print *, 'Electric field file input cadence (s):  ',cfg%dtE0
    print *, 'Electric field file input source directory:  ' // cfg%E0dir
  end if

  if (cfg%flagglow==1) then
    print *, 'GLOW enabled for auroral emission calculations.'
    print *, 'GLOW electron transport calculation cadence (s): ', cfg%dtglow
    print *, 'GLOW auroral emission output cadence (s): ', cfg%dtglowout
  end if
end if

!! default values
lid2in = -1  !< sentinel

do i = 3,argc
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
  case ('-manual_grid')
    call get_command_argument(i+1, argv)
    read(argv,*) lid2in
    call get_command_argument(i+2, argv)
    read(argv,*) lid3in
  end select
end do

end subroutine cli


subroutine help_cli(git_revision)
character(*), intent(in) :: git_revision

print '(/,A,/)', 'MAGCALC ' // git_revision
print '(A)', 'by Matthew Zettergren'
print '(A)', 'build system and software engineering by Michael Hirsch'
print '(A)', 'Compiler vendor: '// get_compiler_vendor()
print '(A)', 'Compiler version: ' // compiler_version()
print '(/,A)', 'must specify input directory and fieldpoint file. Example:'
print '(/,A,/)', 'mpiexec -n 4 build/magcalc.bin test2d_fang test2d_fang/fieldpoint'
print '(A)', '-dryrun option allows quick check of first time step'
stop 'EOF: MAGCALC'
end subroutine help_cli

end module magcalc_cli
