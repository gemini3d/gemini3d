module magcalc_cli

use config, only : read_configfile, gemini_cfg
use pathlib, only : assert_file_exists, assert_directory_exists

implicit none (type, external)
private
public :: cli

contains

subroutine cli(myid, lid, cfg, lid2in, lid3in, debug)

integer, intent(in) :: myid, lid
type(gemini_cfg), intent(out) :: cfg
integer, intent(out) :: lid2in, lid3in
logical, intent(inout) :: debug

integer :: argc, i
character(256) :: argv
character(8) :: date
character(10) :: time

logical :: file_exists

argc = command_argument_count()
if (argc < 2) then
  print '(/,A)', 'MAGCALC: by Matthew Zettergren'
  print '(A,/)', 'build system and software engineering by Michael Hirsch'
  print '(A)', 'must specify input directory and fieldpoint file. Example:'
  print '(/,A,/)', 'mpiexec -n 4 build/magcalc.bin test2d_fang test2d_fang/fieldpoint'
  print '(A)', '-dryrun option allows quick check of first time step'
  stop 'EOF: MAGCALC'
  !! stops with de facto "skip test" return code
endif

call get_command_argument(0, argv)
call date_and_time(date,time)
print '(2A,I6,A3,I6,A)', trim(argv), ' Process:  ', myid,' / ',lid-1, ' at ' // date // 'T' // time


!> READ CONFIG FILE FROM OUTPUT DIRECTORY
call get_command_argument(1,argv)
cfg%outdir = trim(argv)
cfg%infile = cfg%outdir//'/inputs/config.nml'
inquire(file=cfg%infile, exist=file_exists)    !needed to deal with ini vs. nml inputs...
if (.not. file_exists) then
  cfg%infile = cfg%outdir//'/inputs/config.ini'
end if

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
  print '(A,/,A,/,A,/,A)', 'gemini.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

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

end module magcalc_cli
