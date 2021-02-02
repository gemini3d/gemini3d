program test_excess_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use autogrid, only : checker=>test_process_number

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

rx2 = [1,1,2,2,2, 4, 8]
rx3 = [1,2,2,3,4, 3, 6]
N   = [1,2,4,6,8,12,48]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=8, lx3all=12, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,3,2,3,6,6,12]
rx3 = [1,2,2,2,4,4,4,8, 8]
N   = [1,2,4,6,8,12,24,48,96]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=12, lx3all=8, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,2,4, 4, 2, 4, 4,44,  4,  4, 22, 44]
rx3 = [1,2,2,3,2, 3, 9, 6, 9, 2, 27, 54, 18, 18]
N   = [1,2,4,6,8,12,18,24,36,88,108,216,396,792]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=44, lx3all=54, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,3,2, 3, 9, 6, 9, 2, 27, 54, 18, 18]
rx3 = [1,2,2,2,4, 4, 2, 4, 4,44,  4,  4, 22, 44]
N   = [1,2,4,6,8,12,18,24,36,88,108,216,396,792]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=54, lx3all=44, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

if (.not. all_ok) error stop 'process_grid_auto mismatch'

print *, 'OK: auto process grid'

end program
