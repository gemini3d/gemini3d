program main

use filesystem
use, intrinsic :: iso_fortran_env, only: stderr=> error_unit

implicit none

valgrind : block

character(:), allocatable :: in, ref, cwd, out

in = "rel"
cwd = get_cwd()
ref = cwd // filesep() // in

out = absolute(in, '')
if(len_trim(out) == 0) error stop "absolute() has empty output"
if (out /= ref) then
  write(stderr, '(a)') "Mismatch: absolute(" // in //") " // out // " /= " // ref
  error stop
endif

end block valgrind

end program
