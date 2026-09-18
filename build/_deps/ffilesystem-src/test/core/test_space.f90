program test_space_available

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit, int64
use filesystem

implicit none

valgrind : block

character(:), allocatable :: s1
integer(int64) :: i64

s1 = parent(getarg(0))

i64 = space_available(s1)

if (i64 == 0) then
  write(stderr, '(a,i0)') "space_available(" // s1 // ") failed: ", i64
  error stop
end if

print *, "space_available (GB): ", real(i64) / 1024**3

i64 = space_capacity(s1)

if (i64 == 0) then
  write(stderr, '(a,i0)') "space_capacity(" // s1 // ") failed: ", i64
  error stop
end if

print *, "space_capacity (GB): ", real(i64) / 1024**3

end block valgrind

end program
