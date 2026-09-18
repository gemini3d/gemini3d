program test_separator

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use filesystem

implicit none

integer :: c = 0

valgrind : block

character(:), allocatable :: b, s1
character :: sep

sep = pathsep()

if(is_windows() .and. .not. sep == ";") then
  write(stderr,*) "ERROR: pathsep: ", sep
  c = c+1
end if

if(.not. is_windows() .and. .not. sep == ":") then
  write(stderr,*) "ERROR: pathsep: ", sep
  c = c+1
end if

print '(a)', "PASSED: path separator"

b = as_posix("")
if (b /= "") then
  write(stderr,*) "ERROR: as_posix empty: " // b, len_trim(b)
  c = c+1
end if

print '(a)', "PASSED: as_posix('')"

if(is_windows()) then

s1 = "a" // char(92) // "b"
b = as_posix(s1)
if (b /= "a/b") then
  write(stderr,*) "ERROR: as_posix(" // s1 // "): " // b
  c = c+1
end if

end if

!> devnull
b = devnull()
if(is_windows()) then
  if (b /= "nul") then
    write(stderr,*) "ERROR: devnull: " // b
    c = c+1
  end if
else
  if (b /= "/dev/null") then
    write(stderr,*) "ERROR: devnull: " // b
    c = c+1
  end if
end if

print '(a)', "devnull: " // b

end block valgrind

if (c /= 0) error stop "failed separator tests"

print '(a)', "OK: separator"

end program
