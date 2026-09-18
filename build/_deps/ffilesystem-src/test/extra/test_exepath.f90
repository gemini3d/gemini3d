program test_path

use, intrinsic :: iso_fortran_env, only: stderr => error_unit
use, intrinsic :: iso_c_binding, only : C_NULL_CHAR

use filesystem

implicit none

integer :: L, Lt, i
valgrind : block

character(:), allocatable :: path

path = exe_path()

Lt = len_trim(path)

if(Lt == 0) error stop "exe_path() is empty"

L = len(path)

if(Lt /= L) then
  write(stderr, *) "exe_path() has trailing spaces: ", Lt, " /= ", L
  error stop
endif

i = index(path, C_NULL_CHAR)
if(i /= 0) then
  write(stderr,'(a,i0)') "exe_path() " // path // " has embedded null at index ", i
  error stop
endif

if(.not. is_file(path)) error stop "exe_path() is not a file: " // path

print '(a)', path

end block valgrind

end program
