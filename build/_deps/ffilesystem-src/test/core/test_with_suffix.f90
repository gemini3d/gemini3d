program test_with_suffix

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use filesystem

implicit none

valgrind : block

character(:), allocatable :: r

r = with_suffix("", ".h5")
if(r /= ".h5") error stop "with_suffix(,.h5) " // r // " /= .h5"

r = with_suffix("foo.h5", "")
if(r /= "foo") error stop "with_suffix(foo.h5,) " // r // " /= foo"

r = with_suffix(".foo.h5", ".txt")
if(r /= ".foo.txt") error stop "with_suffix(.foo.h5,.txt) " // r // " /= .foo.txt"

r = with_suffix(".h5", "")
if(r /= ".h5") error stop "with_suffix(.h5,) " // r // " /= .h5"

end block valgrind

print '(a)', "OK: with_suffix"

end program
