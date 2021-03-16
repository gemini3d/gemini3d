program gemini3d_run
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use hwloc_ifc, only : get_cpu_count
use runner, only : clean_output
use reader, only: get_simsize3
use autogrid, only : grid_auto, max_mpi
use help, only : help_gemini_run
use config, only : get_compiler_vendor
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, compiler_version

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu, argc
character(1000) :: buf
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra, git_revision
logical :: exists, plan

plan = .false.

git_revision = "@git_rev@"

argc = command_argument_count()

call get_command_argument(1, buf, status=i)
if (i/=0) call help_gemini_run(git_revision)

if(buf(1:1) == '-') then
!! not running sim, checking parameters
  select case(buf)
  case ('-h', '-help')
    call help_gemini_run(git_revision)
  case ('-compiler')
    print '(A)', get_compiler_vendor()
  case ('-compiler_version')
    print '(A)', compiler_version()
  case ('-git')
    print '(A)', git_revision
  case default
    write(stderr,*) "unknown option: ", trim(buf)
    call help_gemini_run(git_revision)
  end select

  stop
endif

Ncpu = get_cpu_count()

@apple_m1_workaround@

print '(A,I0)', 'gemini3d.run: detected CPU count: ', Ncpu

path = trim(buf)

do i = 2, argc
  call get_command_argument(i, buf)

  select case (buf)

  case ('-n')
    call get_command_argument(i+1, buf)
    read(buf, '(I6)') Ncpu
  case ('-gemexe')
    call get_command_argument(i+1, buf)
    gem_exe = trim(buf)
  case ('-mpiexec')
    call get_command_argument(i+1, buf)
    mpiexec = trim(buf)
  case ('-dryrun')
    extra = '-dryrun'
  case ('-plan')
    plan = .true.
  end select
end do

!> setup run
call get_simsize3(path // '/inputs', lx1, lx2all, lx3all)

if(Ncpu > 1) then
  lid = max_mpi(lx2all, lx3all, Ncpu)
else
  lid = 1
endif

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

print '(A,I0,A1,I0,A,I0,A1,I0)', 'MPI partition of lx2, lx3: ', lx2all, ' ',lx3all, &
' is lid2, lid3: ',lid2,' ',lid3
print '(A,I0)', 'MPI images: ', lid

if(plan) stop 'gemini3d.run: plan complete'

!> Find gemini.bin, the main program
if(.not.allocated(gem_exe)) gem_exe = 'gemini.bin'
inquire(file=gem_exe, exist=exists)
if(.not.exists) then
  inquire(file=gem_exe // '.exe', exist=exists)
  if(exists) then
    gem_exe = gem_exe // '.exe'
  endif
endif

inquire(file=gem_exe, exist=exists)
if(.not. exists) error stop "please specify path to gemini.bin with '-gemexe path/to/gemini.bin'"

!> check MPIexec
if(.not.allocated(mpiexec)) then
  call get_environment_variable("MPI_ROOT", buf, status=i)
  if (i==0) mpiexec = trim(buf) // "/bin/mpiexec"
endif
exists = .false.
buf = ""
if(allocated(mpiexec)) then
  inquire(file=mpiexec, exist=exists)
  buf = mpiexec
endif
if(.not.exists) then
  write(stderr,"(A)") "WARNING: MPIexec file not found " // trim(buf)
  write(stderr,"(A)") "If simulation hangs or operates incorrectly, specify -mpiexec option or set MPI_ROOT environment variable."
  mpiexec = "mpiexec"
endif

if(.not.allocated(extra)) extra = ''

!> remove old output files
call clean_output(path)

!> run gemini.bin
if(lid > 1) then
  write(buf, '(A1,A,A,I0,1X,A,1X,A,1X,A)') '"',mpiexec, '" -n ', lid, gem_exe, path, extra
else
  write(buf, '(A,1X,A,1X,A)') gem_exe, path, extra
endif
!! quotes are for mpiexec path with spaces
cmd = trim(buf)
print *, cmd
call execute_command_line(cmd, exitstat=i)

if (i/=0) error stop 'gemini.bin run failure'

end program
