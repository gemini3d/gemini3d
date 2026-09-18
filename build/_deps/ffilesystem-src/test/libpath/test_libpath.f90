program test_binpath

use filesystem

implicit none

valgrind : block

character(:), allocatable :: path

path = lib_path()

if(len_trim(path) == 0) error stop "lib_path() is empty"

if (.not. is_file(path)) error stop path // " is not a file"

print '(a)', path

end block valgrind

end program
