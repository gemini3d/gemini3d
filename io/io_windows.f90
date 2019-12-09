submodule (pathlib) io_windows
!! Windows file functions

use, intrinsic :: iso_c_binding, only: c_char, c_long, c_null_char

implicit none

interface
subroutine fullpath_c(absPath, relPath, maxLength) bind(c, name='_fullpath')
!! char *_fullpath(char *absPath, const char *relPath, size_t maxLength)
!! https://docs.microsoft.com/en-us/cpp/c-runtime-library/reference/fullpath-wfullpath?view=vs-2019
import c_char, c_long
character(kind=c_char), intent(in) :: relPath(*)
character(kind=c_char), intent(out) :: absPath(*)
integer(c_long), intent(in) :: maxLength
end subroutine fullpath_c
end interface

contains

module procedure realpath
!! Windows

integer(c_long), parameter :: N = 260
character(kind=c_char):: c_buf(N)
character(N) :: buf
integer :: i

call fullpath_c(c_buf, path // c_null_char, N)

do i = 1,N
  if (c_buf(i) == c_null_char) exit
  buf(i:i) = c_buf(i)
enddo

realpath = trim(buf(:i-1))
end procedure realpath

end submodule io_windows