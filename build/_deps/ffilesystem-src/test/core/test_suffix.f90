program test_suffix

use filesystem

implicit none

valgrind : block

character(:), allocatable :: s1

if(suffix("") /= "") error stop "suffix empty"

s1 = suffix("suffix_name.a.b")

if (s1 /= ".b") error stop "suffix failed: " // s1

s1 = suffix(s1)
if (s1 /= "") error stop "suffix recursive failed: " // s1

end block valgrind

print '(a)', "OK: Fortran test_suffix"

end program
