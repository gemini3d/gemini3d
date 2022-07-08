module gemini_cli

use gemini3d_config, only : gemini_cfg
use gemini3d_sysinfo, only : get_compiler_vendor
use filesystem, only : expanduser
use mpimod, only : mpibreakdown, mpi_cfg
use exe_frontend, only : help_gemini_bin

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
if (i/=0) then
  call help_gemini_bin()
  ierr = mpibreakdown()
  stop 1
endif

select case (argv)
case ('-h', '-help')
  call help_gemini_bin()
  ierr = mpibreakdown()
  stop
case ('-compiler')
  print '(A)', get_compiler_vendor()
  stop
case ('-git')
  print '(A)', cfg%git_revision
  stop
end select


call get_command_argument(0, argv)
call date_and_time(date,time)
print '(2A,I6,A3,I6,A)', trim(argv), ' Process:  ', mpi_cfg%myid,' / ',mpi_cfg%lid-1, ' at ' // date // 'T' // time


!> READ FILE INPUT
call get_command_argument(1, argv, status=i)
if (i/=0) error stop 'bad command line'
cfg%outdir = expanduser(argv)


!! default values
lid2 = -1  !< sentinel

!> simple but effective command line parsing
do i = 2,argc
  call get_command_argument(i,argv)

  select case (argv)
  case ('-d', '-debug')
    debug = .true.
  case ('-dryrun')
    !! this is a no file output test mode that runs one time step then quits
    !! it helps avoid HPC queuing when a simple setup error exists
    cfg%dryrun = .true.
  case ('-nooutput')
    cfg%nooutput = .true.
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
