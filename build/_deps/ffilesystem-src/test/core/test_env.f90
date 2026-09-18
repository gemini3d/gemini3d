program test_env

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none

call test_exists()
print '(a)', "OK fs: exists"

call test_homedir()
print '(a)', "OK fs: homedir"

if (len_trim(get_tempdir()) == 0) error stop "get_tempdir failed"
print '(a)', "OK: get_tempdir: " // get_tempdir()

call test_username()
print '(a)', "OK fs: username"


contains


subroutine test_exists()

if(exists("")) error stop "empty does not exist"

if(.not. exists(get_cwd())) error stop "exists(get_cwd) failed"

end subroutine


subroutine test_homedir()

character(:), allocatable :: h, p, k, buf

if(is_windows()) then
  k = "USERPROFILE"
else
  k = "HOME"
end if

buf = getenv(k)
if (len_trim(buf) == 0) then
  print '(a)', "env var " // k // " not set"
else
  print '(a)', "getenv: " // k // " = " // trim(buf)
end if

h = get_homedir()
if (len_trim(h) == 0) error stop "get_homedir failed"

p = get_profile_dir()
if (len_trim(p) == 0) error stop "get_profile_dir failed"

print '(a)', "OK: get_homedir: " // h

end subroutine


subroutine test_username()

character(:), allocatable :: u

u = get_username()
if (len_trim(u) == 0) error stop "get_username failed"

print '(a)', "OK: get_username: " // u

end subroutine


end program
