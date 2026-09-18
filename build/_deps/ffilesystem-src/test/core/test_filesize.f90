program test_filesize

use, intrinsic :: iso_fortran_env, only: int64, stderr=>error_unit
use filesystem

implicit none

valgrind : block

integer :: u, d(10)
integer(int64) :: i64

character(:), allocatable :: s1

s1 = parent(getarg(0)) // "/test_filesize.dat"

print '(a)', "file_size path: " // trim(s1)

d = 0

open(newunit=u, file=s1, status="replace", action="write", access="stream")
! writing text made OS-specific newlines that could not be suppressed
write(u) d
close(u)

i64 = file_size(s1)
print '(a, i0)', "filesize (bytes): ", i64
if (i64 /= size(d)*storage_size(d)/8) error stop "file_size() mismatch"
if (i64 /= file_size(s1)) error stop "file_size() mismatch"

end block valgrind

end program
