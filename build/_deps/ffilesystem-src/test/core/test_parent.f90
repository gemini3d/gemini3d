program main

use, intrinsic :: iso_fortran_env, only : stderr => error_unit

use filesystem

implicit none

if (parent("") /= ".") error stop "FAIL: parent()"

if(parent("a/b/c") /= "a/b") error stop "FAIL: parent()"

end program
