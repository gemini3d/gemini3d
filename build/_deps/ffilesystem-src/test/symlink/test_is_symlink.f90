program main

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

use filesystem

implicit none

if(is_symlink("not-exist-file")) error stop "is_symlink() should be false for non-existent file"
if(is_symlink("")) error stop "is_symlink('') should be false"

print '(a)', "OK: is_symlink"

end program
