program test_excess_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use mpimod, only : checker=>test_process_number

implicit none (type, external)

integer, parameter :: N(*) = [1,2,3,4,5,6,8,18,28,64] !< number of fake CPU

!> 2D
call checker(N, lx2all=1, lx3all=40, &
  rx2=[1,1,1,1,1,1,1,1,1,1], &
  rx3=[1,2,3,4,5,6,8,18,28,40])
call checker(N, lx2all=40, lx3all=1, &
  rx2=[1,1,1,1,1,1,1,1,1,1], &
  rx3=[1,2,3,4,5,6,8,18,28,40])

!> 3D
!> NOTE: these are not the most efficient CPU use possible, we
!> can make the auto-grid algorithm more CPU-effective.
call checker(N, lx2all=8, lx3all=12, &
  rx2=[1,1,1,2,2,2,2,4,4,4], &
  rx3=[1,2,3,2,2,3,3,3,3,3])
call checker(N, lx2all=12, lx3all=8, &
  rx2=[1,1,1,2,2,2,2,2,2,2], &
  rx3=[1,2,2,2,2,2,4,4,4,4])

call checker(N, lx2all=6, lx3all=10, &
  rx2=[1,1,1,1,1,1,1,2,2,2], &
  rx3=[1,2,2,2,5,5,5,5,5,5])
call checker(N, lx2all=10, lx3all=6, &
  rx2=[1,1,1,1,1,2,2,2,2,2], &
  rx3=[1,2,3,3,3,3,3,3,3,3])

call checker(N, lx2all=4, lx3all=4, &
  rx2=[1,1,1,2,2,2,2,2,2,2], &
  rx3=[1,2,2,2,2,2,2,2,2,2])

print *, 'OK: auto process grid'

end program
