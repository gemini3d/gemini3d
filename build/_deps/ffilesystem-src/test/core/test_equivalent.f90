program test_canon

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use filesystem

implicit none

call test_same_file()
print *, "OK: same_file"

contains

subroutine test_same_file()

character(:), allocatable :: c1, s1, s2

c1 = getarg(0)

print '(a)', "own command: " // trim(c1)
print '(a)', "current directory: " // get_cwd()

s1 = file_name(c1)
s2 = "./" // s1

if(.not. is_file(s2)) then
  write(stderr, '(a)') "ERROR: file name: " // s2 // " is not a file"
  error stop
end if

if(.not. is_file(s1)) then
  write(stderr, '(a)') "ERROR: file name: " // s1 // " is not a file"
  error stop
end if

print '(a)', "Testing same_file(" // s1 // ", " // s2 // ")"

if(.not. same_file(s1, s2)) then
  write(stderr, '(a)') "ERROR: same_file(" // s1 // " /= " // s2 // ")"
  error stop
end if

s1 = parent(s1)
s2 = "./" // s1

print '(a)', "Testing same_file(" // s1 // ", " // s2 // ")"

if(.not. same_file(s1, s2)) then
  write(stderr, '(a)') "ERROR: same_file(" // s1 // " /= " // s2 // ")"
  error stop
end if

s1 = ".."
s2 = parent(get_cwd())
print '(a)', "Testing same_file(" // s1 // ", " // s2 // ")"
if (.not. same_file(s1, s2)) error stop "ERROR: same_file(" // s1 // " /= " // s2 // ")"

if (.not. same_file(".", "./")) error stop 'ERROR: same_file(., ./)'

!> hardlink resolves to the same file
s1 = "."
s2 = get_cwd()
print '(a)', "Testing same_file(" // s1 // ", " // s2 // ")"
if(.not. same_file(s1, s2)) error stop "ERROR: same_file(" // s1 // " /= " // s2 // ")"

if(same_file("not-exist-same", "not-exist-same")) error stop 'ERROR: same_file(not-exist-same, not-exist-same)'

end subroutine test_same_file

end program
