program test_setenv

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none

valgrind : block
character(:), allocatable :: k, v, buf

k = "TEST_ENV_VAR"
v = "test_value"

call setenv(k, v)
buf = getenv(k)

if (buf /= v) then
  write(stderr, '(a)') "ERROR: setenv/getenv failed: " // buf // " /= " // v
  error stop
end if

end block valgrind

end program
