program test_stem

use filesystem

implicit none


valgrind : block
character(:), allocatable :: s

s = stem("")
if(s /= "") error stop "stem empty: " // s

s = stem("stem.a.b")
if (s /= "stem.a") error stop "%stem failed: " // s
s = stem(s)
if (s /= "stem") error stop "stem nest failed: " // s

end block valgrind

print '(a)', "OK: stem"

end program
