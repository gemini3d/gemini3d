program test

use filesystem

implicit none

! is safe name expects ONLY the filename, not the path

valgrind : block

character(:), allocatable :: s
logical :: ok

s = "test/re/"

if(is_safe_name(s)) error stop s // " is not a safe name"

s = "test/re"
if(is_safe_name(s)) error stop s // " is not a safe name"

s = "hi."
ok = is_safe_name(s)
if(is_windows() .and. ok) error stop "hi. is not a safe name on windows"
if(.not. is_windows() .and. .not. ok) error stop "hi. is a safe name on non-windows"

s = "hi there";
if(is_safe_name(s)) error stop s // " is not a safe name--no spaces allowed"

end block valgrind

print '(a)', "OK: is_safe_name"

end program
