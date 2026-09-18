program main

use filesystem

implicit none

character(*), parameter :: s1 = "a/b/c"

if (.not. is_subdir(s1, "a/b")) error stop "a/b/c is subdir of a/b"
if (.not. is_subdir(s1, "a/b/")) error stop "a/b/c is subdir of a/b/"
if (.not. is_subdir(s1, "a")) error stop "a/b/c is subdir of a"

print '(a)', "OK: is_subdir()"

end program
