program test

use, intrinsic :: iso_fortran_env, only: stderr => error_unit

use filesystem

implicit none

valgrind : block

character(:), allocatable :: fsname

if(is_windows()) then
  fsname = filesystem_type(getenv("SYSTEMDRIVE"))
else
  fsname = filesystem_type("/")
end if

if(len_trim(fsname) == 0) then
  write(stderr, '(a)') "Unknown filesystem type, see type ID in stderr to update fs_get_type()"
  error stop 77
end if

print '(a)', "OK: filesystem_type: " // fsname

end block valgrind

end program
