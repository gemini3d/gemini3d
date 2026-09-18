program test_join

use filesystem

valgrind : block

character(:), allocatable :: p

p = join("", "")
if(p /= "") error stop "join empty: " // p

p = join("/", "")
if(p /= "/") error stop "join(/,): " // p

p = join("", "/")
if(p /= "/") error stop "join('',/): " // p

p = join("ab/cd", "/ef")
if (p /= "/ef") error stop "join(ab/cd, /ef): " // p

end block valgrind

print '(a)', "OK: Fortran join"

end program
