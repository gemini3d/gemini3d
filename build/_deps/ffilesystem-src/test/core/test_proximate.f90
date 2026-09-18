program main

use filesystem

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none

if(proximate_to("", "") /= "") error stop "proximate_to empty failed"

if(proximate_to("Hello", "Hello/") /= ".") error stop "proximate_to same dir failed"

end program
