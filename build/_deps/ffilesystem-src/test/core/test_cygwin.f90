program main

use filesystem
use, intrinsic :: iso_fortran_env, only : stderr => error_unit

implicit none

valgrind : block

character(:), allocatable :: cp, wp, buf

if(.not.is_cygwin()) then
  write(stderr, '(a)') "this test only for Cygwin"
  stop 77
end if

buf = getenv("SYSTEMDRIVE")
if (len_trim(buf) == 0) then
  write(stderr, '(a)') "SYSTEMDRIVE not found"
  stop 77
end if

cp = to_cygpath(buf)
print '(a)', "to_cygpath(" // buf // "): " // cp
if (len_trim(cp) == 0) error stop "to_cygpath(systemdrive) should not be empty"

if (.not. is_dir(cp)) error stop "systemdrive should be an existing directory"

buf = getenv("HOME")
if (len_trim(buf) == 0) then
  write(stderr, '(a)') "HOME not found"
  stop 77
end if

wp = to_winpath(buf)
print '(a)', "to_winpath(" // buf // "): " // wp

end block valgrind

print '(a)', "OK: Cygwin to_cygpath, to_winpath"

end program main
