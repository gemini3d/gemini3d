program main

use filesystem

implicit none

valgrind : block

character(:), allocatable :: p

p = file_name("")
if(p /= "") error stop "FAILED filename('')"

p = file_name("a/b/c")
if(p /= "c") error stop "FAILED filename('a/b/c') " // p

end block valgrind

if (has_filename("")) error stop "FAILED has_filename('')"
if (has_filename("/")) error stop "FAILED has_filename('/')"
if (.not. has_filename(".")) error stop "FAILED has_filename('.')"
if (has_filename("./")) error stop "FAILED has_filename('./')"

print '(a)', "OK: Fortran filename()"

end program
