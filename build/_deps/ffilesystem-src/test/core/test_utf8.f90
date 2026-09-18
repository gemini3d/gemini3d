program utf8_text

use filesystem

implicit none

call test_utf8()

contains

subroutine test_utf8()

character(4), parameter :: smiley = "😀", wink = "😉"
character(6), parameter :: hello = "你好"

character(4) :: u1
character(6) :: u2
character(:), allocatable :: s1

u1 = file_name("./" // smiley)
print '(a)', u1
if (u1 /= smiley) error stop "ERROR: filename UTF8 smiley: " // u1 // " " // smiley

u1 = file_name("./" // wink)
print '(a)', u1
if (u1 /= wink) error stop "ERROR: filename UTF8 wink: " // u1 // " " // wink

u2 = file_name("./" // hello)
print '(a)', u2
if (u2 /= hello) error stop "ERROR: filename UTF8 hello: " // u2 // " " // hello

!> test C allocation for canonical()
s1 = canonical(".")
print '(a)', "canonical(" // "." // ")" // " " // s1


u1 = canonical(smiley)
print '(a)', "canonical(" // smiley // ")" // " " // u1
if (u1 /= smiley) error stop "ERROR: canonical UTF8 smiley: " // u1 // " "// smiley

u1 = canonical(wink)
print '(a)', "canonical(" // wink // ")" // " " // u1
if (u1 /= wink) error stop "ERROR: canonical UTF8 wink: " // u1 // " " // wink

u2 = canonical(hello)
print '(a)', "canonical(" // hello // ")" // " " // u2
if (u2 /= hello) error stop "ERROR: canonical UTF8 hello: " // u2 // " " // hello

deallocate(s1)

end subroutine

end program
