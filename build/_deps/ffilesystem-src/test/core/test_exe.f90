program test_exe

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none


call test_not_exist()
print '(a)', "PASSED: test_not_exist"

call test_is_exe()
print '(a)', "PASSED: test_exist"

contains

subroutine test_not_exist()

character(:), allocatable :: s1

!> empty file
if(is_file("")) error stop "FAILED:test_exe: is_file('') should be false"
if(is_exe("")) error stop "FAILED:test_exe: is_exe('') should be false"

s1 = get_permissions("")
if(len_trim(s1) /= 0) error stop "FAILED:test_exe: get_permissions('') should be empty: " // s1

!> not exist file
s1 = "not-exist-file"
if (is_file(s1)) error stop "FAILED:test_exe: not-exist-file should not exist."
if (is_exe(s1)) error stop "FAILED:test_exe: non-exist-file cannot be executable"
if(len_trim(get_permissions(s1)) /= 0) error stop "FAILED:test_exe: get_permissions('not-exist') should be empty"

end subroutine test_not_exist


subroutine test_is_exe()

character(*), parameter :: noexe = "test_noexe"

character(:), allocatable :: exe

logical :: ok

if(is_wsl() > 0 .and. filesystem_type(".") == "v9fs") then
  print '(a)', "XFAIL:test_exe: WSL with VFS does not support permissions."
  stop 77
end if

exe = getarg(0)
if(len_trim(exe) == 0 ) error stop "could not get own program name"

print '(a)', "test_is_exe: touch(" // noexe // ")"
call touch(noexe)

print '(a)', "set_permissions(" // trim(noexe) // ", executable=.false.)"
call set_permissions(noexe, executable=.false.)

print '(a)', "permissions: " // trim(exe) // " = " // get_permissions(exe)

if (.not. is_exe(exe)) then
  write(stderr,'(a)') "FAILED:test_exe: " // trim(exe) // " is not executable."
  error stop
end if

print '(a)', "permissions: " // trim(noexe) // " = " // get_permissions(noexe)

if (is_exe(noexe)) then
  if(is_windows()) then
    write(stderr,'(a)') "XFAIL:test_exe: " // trim(noexe) // " is executable on Windows."
  else
    write(stderr,'(a)') "FAILED:test_exe: " // trim(noexe) // " is executable and should not be."
    error stop
  end if
end if

print '(a)', "remove(" // trim(noexe) // ")"
call remove(noexe, ok)

end subroutine test_is_exe

end program
