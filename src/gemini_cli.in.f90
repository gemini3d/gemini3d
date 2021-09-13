module gemini_cli

use mpi, only : mpi_init
use config, only : read_configfile, gemini_cfg, get_compiler_vendor
use pathlib, only : expanduser
use mpimod, only : mpisetup, mpibreakdown, mpi_cfg
use exe_frontend, only : help_gemini_bin
use gemini_init, only : find_config, check_input_files

implicit none (type, external)
private
public :: cli

contains


subroutine cli(cfg, lid2, lid3, debug)

type(gemini_cfg), intent(out) :: cfg
integer, intent(out) :: lid2, lid3
logical, intent(inout) :: debug

integer :: argc, i, ierr
character(512) :: argv
character(8) :: date
character(10) :: time

cfg%git_revision = "@git_rev@"

argc = command_argument_count()

call get_command_argument(1, argv, status=i)
if (i/=0) call help_gemini_bin(cfg%git_revision)

select case (argv)
case ('-h', '-help')
  call help_gemini_bin(cfg%git_revision)
case ('-compiler')
  print '(A)', get_compiler_vendor()
  stop
case ('-git')
  print '(A)', cfg%git_revision
  stop
end select

!> INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpi_init(ierr)
if (ierr/=0) error stop 'mpimod: mpi_init'
!! keep mpi_init out of mpisetup for libgemini, where the calling program does MPI_init
call mpisetup()

if(mpi_cfg%lid < 1) error stop 'number of MPI processes must be >= 1. Was MPI initialized properly?'

call get_command_argument(0, argv)
call date_and_time(date,time)
print '(2A,I6,A3,I6,A)', trim(argv), ' Process:  ', mpi_cfg%myid,' / ',mpi_cfg%lid-1, ' at ' // date // 'T' // time


!> READ FILE INPUT
call get_command_argument(1, argv, status=i)
if (i/=0) error stop 'bad command line'
cfg%outdir = expanduser(argv)

call find_config(cfg)

call read_configfile(cfg, verbose=.false.)

call check_input_files(cfg)

!! default values
lid2 = -1  !< sentinel

!> simple but effective command line parsing
do i = 2,argc
  call get_command_argument(i,argv)

  select case (argv)
  case ('-h', '-help')
    ierr = mpibreakdown()
    if (mpi_cfg%myid == 0) call help_gemini_bin(cfg%git_revision)
    stop
  case ('-compiler')
    ierr = mpibreakdown()
    if (mpi_cfg%myid==0) print '(A)', get_compiler_vendor()
    stop
  case ('-git')
    ierr = mpibreakdown()
    if (mpi_cfg%myid==0) print '(A)', cfg%git_revision
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
    if(ierr/=0) error stop '-manual_grid lx2 lx3 parameters are required. lx2 missing'
    read(argv,*) lid2
    call get_command_argument(i+2, argv, status=ierr)
    if(ierr/=0) error stop '-manual_grid lx2 lx3 parameters are required. lx3 missing'
    read(argv,*) lid3
  end select
end do

end subroutine cli

end module gemini_cli
