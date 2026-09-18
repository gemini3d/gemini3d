program test_normal

use filesystem

implicit none

valgrind : block

character(:), allocatable :: p

p = normal("")
if(p /= ".") error stop "normal() " // p

p = normal("a/b/.")
if(p /= "a/b") error stop "normal(a/b/.) " // p

end block valgrind

print '(a)', "OK: normal"

end program
