program test_touch

use, intrinsic :: iso_c_binding, only : C_LONG_LONG
use filesystem

implicit none

integer(C_LONG_LONG) :: t0, t1


! call touch("")  !< error stops

valgrind : block

character(:), allocatable :: fn, dir

dir = parent(getarg(0))

!> touch absolute path
fn = dir // "/test_fileop.h5"
print '(a)', "touch(" // fn // ")"
call touch(fn)
call assert_is_file(fn)
t0 = get_modtime(fn)
if(t0 < 0) error stop "absolute path: Could not get modtime - is negative integer"
if(t0 == 0) error stop "absolute Path: Could not get modtime - integer 0 value"

print '(a,i0)', "PASSED: touch absolute path: " // fn // " mod time (posix): ", t0


!> touch relative path
print '(a)', "touch relative path, cwd: " // get_cwd()

fn = "test_fileop.h5"

call touch(fn)
call assert_is_file(fn)
t0 = get_modtime(fn)

if(t0 < 0) error stop "relative path: Could not get modtime - is negative integer"
if(t0 == 0) error stop "relative path: Could not get modtime - integer 0 value"

if(.not. set_modtime(fn)) error stop "Could not set modtime"
t1 = get_modtime(fn)

if(t1 < t0) error stop "get_modtime: Could not set modtime"

if(set_modtime('not-exist-tile')) error stop "set_modtime: Could set modtime for non-existing file"

print '(a)', "PASSED: touch relative path: " // fn

end block valgrind

print '(a,1x,i0)', "OK: touch. modtime: ", t0

end program
