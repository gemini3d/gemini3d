program main

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none

if(is_wsl() > 0 .and. filesystem_type(".") == "v9fs") then
  write(stderr, '(a)') "ERROR: WSL with v9fs does not support file permissions"
  stop 77
end if

call test_set_permissions()
print '(a)', "OK: set_permission"

call test_get_permissions()
print '(a)', "OK: get_permission"


contains

subroutine test_get_permissions()

character(9) :: p

character(*), parameter :: reada="readable.txt", noread="not-readable.txt", nowrite="not-writable.txt"
logical :: ok

p = get_permissions("")
if(len_trim(p) /= 0) then
    write(stderr, '(a)') "get_permissions('') should be empty: " // p
    error stop
end if

!> readable
call touch(reada, ok)
if (.not. ok) then
    write(stderr, '(a)') "ERROR: could not create file " // reada
    error stop 77
end if

call set_permissions(reada, readable=.true.)

p = get_permissions(reada)
print '(a)', "Permissions for " // trim(reada)// ": "// p

if (len_trim(p) == 0) error stop "get_permissions('"//trim(reada)//"') should not be empty"

if (p(1:1) /= "r") then
    write(stderr, '(a)') "ERROR: "//trim(reada)//" should be readable"
    error stop
end if
if(.not. is_readable(reada)) error stop trim(reada)//" should be readable"

if(.not. exists(reada)) error stop trim(reada)//" should exist"

if(.not. is_file(reada)) error stop trim(reada)//" should be a file"

!! for Ffilesystem, even non-readable files "exist" and are "is_file"

call touch(noread, ok)
if(.not. ok) then
    write(stderr, '(a)') "ERROR: could not create file " // noread
    error stop 77
end if

call set_permissions(noread, readable=.false.)

p = get_permissions(noread);
print '(a)', "Permissions for " // trim(noread)// ": "// p

if (len_trim(p) == 0) error stop "get_permissions('"//trim(noread)//"') should not be empty"

if(.not. (is_windows() .or. is_cygwin())) then
if (p(1:1) == "r") error stop trim(noread)//" should not be readable " // p
if (index(p, "r") == 0) then
  if (is_readable(noread)) error stop "is_readable: "//trim(noread)//" should not be readable"
end if
end if

if (.not. exists(noread)) error stop trim(noread)//" should exist"
if (.not. is_file(noread)) error stop trim(noread)//" should be a file"

!> writable
if(.not. is_file(nowrite)) call touch(nowrite)
call set_permissions(nowrite, writable=.false.)

if (.not. is_writable(".")) error stop ". should be writable"

p = get_permissions(nowrite);
print '(a)', "Permissions for " // trim(nowrite)// ": "// p
if (len_trim(p) == 0) error stop "get_permissions('"//trim(nowrite)//"') should not be empty"

if (p(2:2) == "w") then
    write(stderr, '(a)') "get_permissions: " //trim(nowrite)//" should not be writable"
    if(.not. (is_windows() .or. is_cygwin())) error stop
end if

if(index(p, "w") == 0 .and. .not. is_admin() .and. is_writable(nowrite)) &
  error stop "is_writable: " // trim(nowrite)//" should not be writable"

if (.not. exists(nowrite)) error stop trim(nowrite)//" should exist"

if (.not. is_file(nowrite)) error stop trim(nowrite)//" should be a file"

call remove(reada, ok)
call remove(noread, ok)
call remove(nowrite, ok)

end subroutine test_get_permissions


subroutine test_set_permissions()

character(9) :: perm
logical :: ok

character(:), allocatable :: exe, noexe

exe = join(get_tempdir(), "yes_exe")
noexe = join(get_tempdir(), "no_exe")

!> chmod(.true.)

call touch(exe)
if(.not. is_file(exe)) error stop "ERROR: not a file: " // exe

perm = get_permissions(exe)
print '(a)', "permissions before set_permissions(exe=true) " // perm

call set_permissions(exe, executable=.true.)

call set_permissions(exe, executable=.false., ok=ok)
if (.not. ok) error stop "ERROR: set_permissions(exe=.false.) failed"

call set_permissions(exe, executable=.true., ok=ok)

perm = get_permissions(exe)
print '(a)', "permissions after set_permissions(exe=true): " // perm

if (.not. is_exe(exe)) then
    write(stderr,'(a)') "ERROR: is_exe() did not detect executable file " // trim(exe)
    if(.not. (is_windows() .or. is_cygwin())) error stop
end if

if(.not. (is_windows() .or. is_cygwin())) then
if(perm(3:3) /= "x") then
    write(stderr,'(a)') "ERROR: get_permissions() " // trim(exe) // " is not executable"
    error stop
end if
end if

!> chmod(.false.)

call touch(noexe)
if(.not. is_file(noexe)) then
    write(stderr,'(a)') "ERROR: " // trim(noexe) // " is not a file."
    error stop
end if

perm = get_permissions(noexe)
print '(a)', "permissions: " // trim(noexe) // " = " // perm

call set_permissions(noexe, executable=.false., ok=ok)
if (.not. ok) error stop "ERROR: set_permissions(exe=.false.) failed"

if(.not. (is_windows() .or. is_cygwin())) then
!~ Windows file system is always executable to stdlib.

    if (is_exe(noexe)) error stop "ERROR:t did not detect non-executable file."

    if(perm(3:3) /= "-") then
    write(stderr,'(a)') "ERROR:get_permissions() " // trim(noexe) // " is executable"
    error stop
    end if

end if

end subroutine test_set_permissions

end program
