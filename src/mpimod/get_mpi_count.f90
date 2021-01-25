program get_mpi_count
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use autogrid, only : grid_auto, max_mpi
use reader, only: get_simsize3

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu, argc
character(1000) :: buf
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra
logical :: exists

argc = command_argument_count()

if(argc < 2) error stop 'specify top-level simulation path and cpu count'

call get_command_argument(1, buf, status=i)
path = trim(buf)

call get_command_argument(2, buf)
read(buf, '(I6)') Ncpu

call get_simsize3(path // '/inputs', lx1, lx2all, lx3all)

lid = max_mpi(lx2all, lx3all, Ncpu)

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

if(argc < 3) then
  write(buf, '(I6,A)') lid, ' CPU cores would be used. Give path to gemini.bin to run simulation'
  stop trim(buf)
end if

call get_command_argument(3, buf)
gem_exe = trim(buf)
inquire(file=gem_exe, exist=exists)
if(.not.exists) error stop gem_exe // ' is not a file.'

if (argc >= 4) then
  call get_command_argument(4, buf)
  mpiexec = trim(buf)
  inquire(file=mpiexec, exist=exists)
  if(.not.exists) error stop mpiexec // ' is not a file.'
else
  mpiexec = 'mpiexec'
endif

if (argc >= 5) then
  call get_command_argument(5, buf)
  extra = trim(buf)
else
  extra = ''
endif

!> run gemini.bin
write(buf, '(A1,A,A,I6,1X,A,1X,A,1X,A)') '"',mpiexec, '" -n ', lid, gem_exe, path, extra
cmd = trim(buf)
print *, cmd
call execute_command_line(cmd, exitstat=i)

if(i/=0) error stop 'gemini.bin run failure'

end program
