program test_realpath_cfi

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use, intrinsic :: iso_c_binding, only : C_INT
use filesystem, only : max_path

implicit none

interface
  integer(C_INT) function realpath_cfi(path, result) bind(C, name="realpath_cfi")
    import C_INT
    character(*), intent(in) :: path
    character(*), intent(out) :: result
  end function realpath_cfi
end interface

character(:), allocatable :: out
character(4) :: tiny
integer(C_INT) :: rc

allocate(character(max_path()) :: out)

rc = realpath_cfi(".", out)
if (rc <= 0) then
  write(stderr, '(a,1x,i0)') "FAILED: realpath_cfi('.') expected > 0 rc, got", rc
  error stop
end if

rc = realpath_cfi(".   ", out)
if (rc <= 0) then
  write(stderr, '(a,1x,i0)') "FAILED: realpath_cfi('.   ') expected > 0 rc, got", rc
  error stop
end if

rc = realpath_cfi("not-exist-realpath-cfi/a/b/c", out)
if (rc /= 0) then
  write(stderr, '(a,1x,i0)') "FAILED: realpath_cfi(nonexistent) expected rc=0, got", rc
  error stop
end if

rc = realpath_cfi(".", tiny)
if (rc /= -5_C_INT) then
  write(stderr, '(a,1x,i0)') "FAILED: realpath_cfi(truncation) expected rc=-5, got", rc
  error stop
end if

print *, "OK: realpath_cfi contract"

end program test_realpath_cfi