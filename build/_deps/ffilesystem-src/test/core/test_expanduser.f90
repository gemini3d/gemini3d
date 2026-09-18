program test_expanduser

use, intrinsic:: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none


if(expanduser("") /= "") error stop "expanduser blank failed"
if(expanduser(".") /= ".") error stop "expanduser dot failed: " // expanduser(".")

if(expanduser("~P") /= "~P") error stop "expanduser ~P failed: " // expanduser("~P")

valgrind: block

character(:), allocatable :: s1, s2, s3

!> does expanduser() get homedir correctly
s1 = expanduser("~")
if(s1 /= get_homedir()) then
   write(stderr, '(a)') "expanduser(~) failed: " // s1 // " /= " // get_homedir()
   error stop
end if

!> idempotent
s2 = expanduser(s1)
if (s1 /= s2) then
  write(stderr, '(a)') "expanduser() idempotent failed: " // s1 // " /= " // s2
  error stop
end if

!> separator
s1 = get_homedir() // "/.."

!> double dot
s3 = "~/.."
s2 = expanduser(s3)
if (s2 /= s1) error stop "expanduser(" // s3 // ") failed: " // s2 // " /= " // s1
print *, "PASS: expanduser("//s3//")  ", s2


end block valgrind

print *, "OK: filesystem: expanduser  ", expanduser("~")

end program
