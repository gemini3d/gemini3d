submodule (pathlib) posix_io
!! for systems supporting POSIX file functions
!! Linux, not MacOS or MinGW

use, intrinsic :: iso_c_binding, only: c_char, c_null_char

implicit none

interface
subroutine realpath_c(path, resolved_path) bind(c, name='realpath')
!! char *realpath(const char *path, char *resolved_path)
!! http://man7.org/linux/man-pages/man3/realpath.3.html
import c_char
character(kind=c_char), intent(in) :: path(*)
character(kind=c_char), intent(out) :: resolved_path(*)
end subroutine realpath_c
end interface

contains

module procedure realpath
!! for POSIX systems (not Windows)
!!
!! N=4096 may be a maximum, can be shorter as per expected maximum path length.
!!
!! * Fortran -> character scalar to C as null-terminated array with pointer to first element
!! * <- C returns a single-character array with a null terminator as the last element.
!!
!! We use the do loop to copy from the character array into a scalar character, stopping at the null terminator
!! this is probably at least as effective as using findloc().

integer, parameter :: N = 2048
character(kind=c_char):: c_buf(N)
character(N) :: buf
integer :: i

call realpath_c(path // c_null_char, c_buf)

do i = 1,N
  if (c_buf(i) == c_null_char) exit
  buf(i:i) = c_buf(i)
enddo

realpath = trim(buf(:i-1))

end procedure realpath

end submodule posix_io