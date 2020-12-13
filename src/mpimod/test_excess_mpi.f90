program test_excess_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use mpimod, only : checker=>test_process_number

implicit none (type, external)

integer :: i
logical :: ok, all_ok
integer, allocatable, dimension(:) :: rx2, rx3, N

all_ok = .true.

!> 2D
rx2 = [1,1,1,1,1, 1, 1, 1]
rx3 = [1,2,4,5,8,10,20,40]
N   = [1,2,4,5,8,10,20,40]
do i = 1,size(N)
  ok = checker(N(i), lx2all=1, lx3all=40, rx2=rx2(i), rx3=rx3(i))
  if (.not. ok) all_ok = .false.
  ok = checker(N(i), lx2all=40, lx3all=1, rx2=rx2(i), rx3=rx3(i))
  if (.not. ok) all_ok = .false.
enddo

!> 3D

rx2 = [1,2,4,2,8, 4, 8]
rx3 = [1,1,1,3,1, 3, 6]
N   = [1,2,4,6,8,12,48]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=8, lx3all=12, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,2,4,6,4,12,12,12,12]
rx3 = [1,1,1,1,2,1 , 2, 4, 8]
N   = [1,2,4,6,8,12,24,48,96]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=12, lx3all=8, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,2,4,2,4, 4, 2, 4, 1, 4,44,  4,  4, 44, 44]
rx3 = [1,1,1,3,2, 3, 9, 6,27, 9, 2, 27, 54,  9, 18]
N   = [1,2,4,6,8,12,18,24,27,36,88,108,216,396,792]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=44, lx3all=54, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,2,2,6,2, 6,18, 6,27,18, 2, 54, 54, 18, 18]
rx3 = [1,1,2,1,4, 2, 1, 4, 1, 2,44,  2,  4, 22, 44]
N   = [1,2,4,6,8,12,18,24,27,36,88,108,216,396,792]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=54, lx3all=44, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

if (all_ok) stop 'OK: auto process grid'

error stop

end program
