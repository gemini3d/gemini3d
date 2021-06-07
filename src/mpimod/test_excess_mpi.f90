program test_excess_mpi

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use autogrid, only : checker=>test_process_number, max_mpi

implicit none (type, external)

call check_autogrid()

call check_auto_mpi()

contains

subroutine check_autogrid()

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
  if (.not. checker(N(i), lx2all=40, lx3all=1, rx2=rx3(i), rx3=rx2(i))) all_ok = .false.
enddo

!> 3D

rx2 = [1,1,2,2,2, 4, 4]
rx3 = [1,2,2,3,4, 3, 6]
N   = [1,2,4,6,8,12,24]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=8, lx3all=12, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
enddo

rx2 = [1,1,2,3,2, 3]
rx3 = [1,2,2,2,4, 4]
N   = [1,2,4,6,8,12]
do i = 1,size(N)
  if (.not. checker(N(i), lx2all=6, lx3all=8, rx2=rx2(i), rx3=rx3(i))) all_ok = .false.
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

end subroutine check_autogrid


subroutine check_auto_mpi()

integer :: i, j, O
integer,parameter :: &
L(3) = [4,16,40], &  !< length of dimension
N(3) = [4,28,96], & !< number of CPU cores
M(3,3) = reshape(&    !< expected auto-selected number of MPI images
[2, 2, 2, &
 4, 8, 8, &
 4, 20, 20], shape=shape(M), order=[2,1])

if (max_mpi(100,52, 1) /= 1) error stop "gcd(N,1) /= 1"

if (max_mpi(6, 8, 32) /= 12) error stop "max_mpi(6,8,32) /= 12"

do i = 1,size(L)
  do j = 1,size(N)
    O = max_mpi(L(i), 1, N(j))
    if (O /= M(i,j)) then
      write(stderr,*) "ERROR:max_mpi",L(i),N(j),M(i,j),O
      error stop
    endif
    if (M(i,j) /= max_mpi(1, L(i), N(j))) error stop 'max_mpi value error'
  enddo
enddo

print *, "OK: max_mpi"

end subroutine check_auto_mpi

end program
