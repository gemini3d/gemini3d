program gemini3d_run
!! for use from terminal/CMake, computes optimal MPI count
!! for a particular simulation

use exe_frontend, only : cli_parser, get_Ncpu, clean_output
use reader, only: get_simsize3
use autogrid, only : grid_auto, max_mpi

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, lid, lid2, lid3, Ncpu
character(:), allocatable :: path, gem_exe, cmd, mpiexec, extra
logical :: plan
character(1000) :: buf

call cli_parser(plan, Ncpu, path, gem_exe, mpiexec, extra)

Ncpu = get_Ncpu(Ncpu)

!> setup run
call get_simsize3(path // '/inputs/simsize.h5', lx1, lx2all, lx3all)

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
