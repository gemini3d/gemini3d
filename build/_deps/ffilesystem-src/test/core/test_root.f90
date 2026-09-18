program test_root

use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use filesystem

implicit none

if(root("") /= "") error stop "root empty"

if(is_windows()) then
  if(root_name("C:/a/b") /= "C:") error stop "root name should be C:"
else
  if(root("/a/b") /= "/") error stop "relative root should be empty"
  if(root_name("/a/b") /= "") error stop "root name should be empty"
endif

print '(A)', "OK: root: C++"

end program
