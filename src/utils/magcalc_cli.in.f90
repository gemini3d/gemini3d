module magcalc_cli

use, intrinsic :: iso_fortran_env, only : compiler_version
use gemini3d_config, only : read_configfile, gemini_cfg
use gemini3d_sysinfo, only : get_compiler_vendor
use filesystem, only : assert_is_file, assert_is_dir, expanduser
use mpimod, only : mpisetup, mpibreakdown, mpi_cfg
use phys_consts, only : wp
use timeutils, only : dateinc
use exe_frontend, only : help_magcalc_bin

implicit none (type, external)
private
public :: cli

contains

subroutine cli(cfg, lid2, lid3, debug, ymdstart,UTsecstart,ymdend,UTsecend)

type(gemini_cfg), intent(inout) :: cfg
integer, intent(out) :: lid2, lid3
logical, intent(inout) :: debug
integer, dimension(3), intent(out) :: ymdstart
real(wp), intent(out) :: UTsecstart
integer, dimension(3), intent(out) :: ymdend
real(wp), intent(out) :: UTsecend

integer :: argc, i, iarg, ierr
character(256) :: argv
character(8) :: date
character(10) :: time


cfg%git_revision = "@git_rev@"

argc = command_argument_count()

call get_command_argument(1, argv, status=i)
if (i/=0) call help_magcalc_bin()

select case (argv)
case ('-h', '-help')
  call help_magcalc_bin()
case ('-compiler')
  print '(A)', get_compiler_vendor()
  stop
case ('-git')
  print '(A)', cfg%git_revision
  stop
end select

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()

if(mpi_cfg%lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'

call get_command_argument(0, argv)
call date_and_time(date,time)
print '(2A,I6,A3,I6,A)', trim(argv), ' Process:  ', mpi_cfg%myid,' / ', mpi_cfg%lid-1, ' at ' // date // 'T' // time


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
iarg = 2 !< start of -options
if (argc > 1) then
  call get_command_argument(2, argv, status=i)
  if (i/=0) error stop 'please specify fieldpoint file'
  if (argv(1:1) /= "-") then
    cfg%fieldpointfile = trim(argv)
    iarg = 3
  endif
endif

if (iarg == 2) then
  find_pts : block
  logical :: exists
  character(*), parameter :: exts(2) = [character(3) :: "h5", "dat"]
  character(:), allocatable :: loc
  do i = 1,size(exts)
    loc = trim(cfg%outdir // "/inputs/magfieldpoints." // exts(i))
    inquire(file=loc, exist=exists)
    if (exists) then
      cfg%fieldpointfile = loc
      exit find_pts
    endif
  end do

  if (cfg%outdir(1:1) == "-") then
    error stop 'magcalc.bin: not a known CLI option: ' // cfg%outdir
  else
    error stop 'magcalc.bin: could not find magfieldpoints file in ' // cfg%outdir
  endif

  end block find_pts
endif

!! this file contains the field points at which we are computing magnetic perturbations, it will be copied into the output directory

if (mpi_cfg%myid==0) then
  print *, 'Simulation data directory:  ', cfg%outdir
  print *, 'Input config file:  ',cfg%infile
  print *, 'fieldpoint file: ', cfg%fieldpointfile
end if

call read_configfile(cfg)


!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (mpi_cfg%myid==0) then
  call assert_is_file(cfg%indatsize)
  call assert_is_file(cfg%indatgrid)
  call assert_is_file(cfg%indatfile)

  print '(A,I6,A1,I0.2,A1,I0.2)', cfg%infile // ' start year-month-day:  ',cfg%ymd0(1),'-',cfg%ymd0(2),'-',cfg%ymd0(3)
  print '(A51,F10.3)', 'start time:  ',cfg%UTsec0
  print '(A51,F10.3)', 'duration:  ',cfg%tdur
  print '(A51,F10.3)', 'output every:  ',cfg%dtout
  print '(A,/,A,/,A,/,A)', 'magcalc.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

  if(cfg%flagdneu==1) then
    call assert_is_dir(cfg%sourcedir)
    print *, 'Neutral disturbance mlat,mlon:  ',cfg%sourcemlat,cfg%sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',cfg%dtneu
    print *, 'Neutral grid resolution (m):  ',cfg%drhon,cfg%dzn
    print *, 'Neutral disturbance data files located in directory:  ',cfg%sourcedir
  end if

  if (cfg%flagprecfile==1) then
    call assert_is_dir(cfg%precdir)
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',cfg%dtprec
    print *, 'Precipitation file input source directory:  ' // cfg%precdir
  end if

  if(cfg%flagE0file==1) then
    call assert_is_dir(cfg%E0dir)
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
lid2 = -1  !< sentinel

ymdstart = [0,0,0]
UTsecstart = 0
ymdend = [0,0,0]
UTsecend = 0

do i = iarg,argc
  call get_command_argument(i,argv)

  select case (argv)
  case ('-d', '-debug')
    debug = .true.
  case ('-dryrun')
    !! this is a no file output test mode that runs one time step then quits
    !! it helps avoid HPC queuing when a simple setup error exists
    cfg%dryrun = .true.
  case ('-manual_grid')
    call get_command_argument(i+1, argv, status=ierr)
    if(ierr/=0) error stop '-manual_grid lx2 lx3 parameters are required. lx2 missing'
    read(argv, '(I6)') lid2
    call get_command_argument(i+2, argv, status=ierr)
    if(ierr/=0) error stop '-manual_grid lx2 lx3 parameters are required. lx3 missing'
    read(argv, '(I6)') lid3
  case ('-start_time')
    call get_command_argument(i+1, argv, status=ierr)
    if(ierr/=0) error stop '-start_time year month day UTsec parameters are required. year missing'
    read(argv, '(I4)') ymdstart(1)
    call get_command_argument(i+2, argv, status=ierr)
    if(ierr/=0) error stop '-start_time year month day UTsec parameters are required. month missing'
    read(argv, '(I2)') ymdstart(2)
    call get_command_argument(i+3, argv, status=ierr)
    if(ierr/=0) error stop '-start_time year month day UTsec parameters are required. day missing'
    read(argv, '(I2)') ymdstart(3)
    call get_command_argument(i+4,argv, status=ierr)
    if(ierr/=0) error stop '-start_time year month day UTsec parameters are required. UTsec [0..86400) missing'
    read(argv, '(F9.3)') UTsecstart
  case ('-end_time')
    call get_command_argument(i+1, argv, status=ierr)
    if(ierr/=0) error stop '-end_time year month day UTsec parameters are required. year missing'
    read(argv, '(I4)') ymdend(1)
    call get_command_argument(i+2, argv, status=ierr)
    if(ierr/=0) error stop '-end_time year month day UTsec parameters are required. month missing'
    read(argv, '(I2)') ymdend(2)
    call get_command_argument(i+3, argv, status=ierr)
    if(ierr/=0) error stop '-end_time year month day UTsec parameters are required. day missing'
    read(argv, '(I2)') ymdend(3)
    call get_command_argument(i+4, argv, status=ierr)
    if(ierr/=0) error stop '-end_time year month day UTsec parameters are required. UTsec [0..86400) missing'
    read(argv, '(F9.3)') UTsecend
  end select
end do

! have root report start and end time for user
if (mpi_cfg%myid==0) then
  print*, 'Start time requested for magcalc:  ',ymdstart,UTsecstart
  print*, 'End time requested for magcalc:  ',ymdend,UTsecend
end if

end subroutine cli

end module magcalc_cli
