program gemini_runner
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use runner, only : clean_output, get_cpu_count
use reader, only: get_simsize3
use autogrid, only : grid_auto, max_mpi
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu, argc
character(1000) :: buf
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra
logical :: exists

Ncpu = 0

argc = command_argument_count()

if(argc < 1) error stop 'specify top-level simulation path'

call get_command_argument(1, buf)
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
    inquire(file=mpiexec, exist=exists)
    if(.not.exists) error stop mpiexec // ' is not a file.'
  case ('-dryrun')
    extra = '-dryrun'
  end select
end do

if(.not.allocated(gem_exe)) then
  gem_exe = 'gemini.bin'
  inquire(file=gem_exe, exist=exists)
  if(.not.exists) gem_exe = "gemini.bin.exe"
  inquire(file=gem_exe, exist=exists)
  if(.not.exists) error stop "please specify path to gemini.bin with '-gemexe path/to/gemini.bin'"
else
  inquire(file=gem_exe, exist=exists)
  if(.not.exists) error stop gem_exe // ' is not a file.'
endif

if(.not.allocated(mpiexec)) mpiexec = 'mpiexec'

if(.not.allocated(extra)) extra = ''

!> remove old output files
call clean_output(path)

!> setup run
call get_simsize3(path // '/inputs', lx1, lx2all, lx3all)

if(Ncpu == 0) then
  Ncpu = get_cpu_count()
  print '(A,I0)', 'Detected CPU count: ', Ncpu
endif

if(Ncpu > 1) then
  lid = max_mpi(lx2all, lx3all, Ncpu)
else
  lid = 1
endif

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

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
