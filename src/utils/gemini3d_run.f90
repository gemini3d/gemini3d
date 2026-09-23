program gemini3d_run
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use exe_frontend, only : cli_parser, get_Ncpu, quote_argument
use reader, only: get_simsize3
use autogrid, only : grid_auto, max_mpi

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra
logical :: plan
character(20) :: count

call cli_parser(plan, path, gem_exe, mpiexec, extra)

Ncpu = get_Ncpu()

!> setup run
call get_simsize3(path // '/inputs/simsize.h5', lx1, lx2all, lx3all)

lid = max_mpi(lx2all, lx3all, Ncpu)

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

!> JSON format
print '(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)', '{ "lx2": ', lx2all, ', "lx3": ', lx3all, ', "lid2": ', lid2, ', "lid3": ', lid3, &
', "lid": ', lid, ', "Ncpu": ', Ncpu, ' }'

if(.not. plan) then
! we didn't use STOP because that prints "STOP" with Gfortran.

!> run gemini.bin
if(lid > 1) then
  write(count, '(I0)') lid
  cmd = quote_argument(mpiexec) // ' -n ' // trim(count) // ' '
else
  cmd = ''
endif
cmd = cmd // quote_argument(gem_exe) // ' ' // quote_argument(path)
if(len_trim(extra)>0) cmd = trim(cmd) // ' ' // trim(adjustl(extra))
print *, cmd
call execute_command_line(cmd, exitstat=i)

if (i/=0) error stop 'gemini.bin run failure'

endif

end program
