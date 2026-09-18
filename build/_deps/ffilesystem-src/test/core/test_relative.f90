program main

use filesystem

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none

if(relative_to("", "") /= "") error stop "FAIL: relative_to()"

if(relative_to("a/b/c", "a/b") /= "..") error stop "FAIL: relative_to()"

print '(a)', "PASS: relative_to()"

end program
