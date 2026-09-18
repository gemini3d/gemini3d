program test_windows

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none

character(:), allocatable :: s1, s2, s3

s1 = getenv("PROGRAMFILES")
if (len(s1) == 0) then
    write(stderr, '(a)') "Error getting PROGRAMFILES"
    error stop 77
end if

print '(a)', "PROGRAMFILES: " // trim(s1)

s2 = shortname(s1)
print '(a)', trim(s1) // " => " // trim(s2)
if(len_trim(s2) == 0) then
    write(stderr, '(a)') "Error converting long path to short path: " // trim(s2)
    error stop
end if

s3 = longname(s2)
print '(a)', trim(s2) // " => " // trim(s3)
if(len_trim(s3) == 0) then
    write(stderr, '(a)') "Error converting short path to long path: " // trim(s3)
    error stop
end if

if(s1 /= s3) then
    write(stderr, '(a)') "Error: shortname(longname(x)) != x"
    error stop
end if

end program
