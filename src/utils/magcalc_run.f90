program magcalc_run
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use exe_frontend, only : cli_parser, get_Ncpu, quote_argument
use reader, only: get_simsize3
use autogrid, only : grid_auto, max_mpi

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu
character(:), allocatable :: path, exe, cmd, mpiexec, numproc_flag, extra
logical :: plan
character(20) :: count
call cli_parser(plan, path, exe, mpiexec, numproc_flag, extra)

Ncpu = get_Ncpu()

if (Ncpu <= 1) error stop 'Ncpu must be > 1. use mpiexec with magcalc.bin'

!> setup run
call get_simsize3(path // '/inputs/simsize.h5', lx1, lx2all, lx3all)

lid = max_mpi(lx2all, lx3all, Ncpu)

!> checks consistency
call grid_auto(lx2all, lx3all, lid, lid2, lid3)

print '(A,I0,A1,I0,A,I0,A1,I0)', 'MPI partition of lx2, lx3: ', lx2all, ' ',lx3all, &
' is lid2, lid3: ',lid2,' ',lid3
print '(A,I0)', 'MPI images: ', lid

if(plan) stop 'magcalc.run: plan complete'

!> run magcalc.bin
if(lid > 1) then
  write(count, '(I0)') lid
  cmd = quote_argument(mpiexec) // ' ' // trim(numproc_flag) // ' ' // trim(count) // ' '
else
  cmd = ''
endif
cmd = cmd // quote_argument(exe) // ' ' // quote_argument(path)
if(len_trim(extra)>0) cmd = trim(cmd) // ' ' // trim(adjustl(extra))
print *, cmd
call execute_command_line(cmd, exitstat=i)

if (i/=0) error stop 'magcalc.bin run failure'

end program
