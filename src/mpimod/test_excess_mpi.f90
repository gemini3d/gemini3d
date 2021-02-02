program test_excess_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use autogrid, only : checker=>test_process_number

implicit none (type, external)

integer :: i
logical :: all_ok
integer, allocatable, dimension(:) :: rx2, rx3, N

all_ok = .true.

!> 2D
rx2 = [1,1,1,1,1, 1, 1, 1]
rx3 = [1,2,4,5,8,10,20,40]
N   = [1,2,4,5,8,10,20,40]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=1, lx3all=40, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
  if (.not. checker(N(i), lx2all=40, lx3all=1, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

!> 3D

rx2 = [1,1,2,2,2, 4, 4]
rx3 = [1,2,2,3,4, 3, 6]
N   = [1,2,4,6,8,12,24]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=8, lx3all=12, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,5,2, 4, 5, 5, 10]
rx3 = [1,2,2,1,4, 4, 4, 8,  8]
N   = [1,2,4,5,8,16,20,40, 80]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=20, lx3all=16, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,2,4, 4, 2, 4, 4,11,  4, 11, 11, 22]
rx3 = [1,2,2,3,2, 3, 9, 6, 9, 9, 27, 18, 27, 27]
N   = [1,2,4,6,8,12,18,24,36,99,108,198,297,594]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=44, lx3all=54, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,3,2, 3, 9, 6, 9, 9, 27, 18, 27, 27]
rx3 = [1,2,2,2,4, 4, 2, 4, 4,11,  4, 11, 11, 22]
N   = [1,2,4,6,8,12,18,24,36,99,108,198,297,594]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=54, lx3all=44, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

if (.not. all_ok) error stop 'process_grid_auto mismatch'

print *, 'OK: auto process grid'

end program
