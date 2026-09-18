program test

use filesystem

implicit none

valgrind : block

character(:), allocatable :: pwd, p1
logical :: ok

call mkdir("", ok)
if(ok) error stop "mkdir: empty path"

pwd = get_cwd()
if(len_trim(pwd) == 0) error stop "get_cwd: " // pwd

print '(a)', "test_mkdir: pwd " // pwd

call mkdir(pwd)
!! test that existing directory does not cause error

p1 = pwd // "/test-filesystem-dir1"
call mkdir(p1)
if(.not. is_dir(p1)) error stop "mkdir: full-path: " // p1
print '(a)', "PASS: mkdir full-path"

call mkdir(p1, ok=ok)
if(.not. ok) error stop "mkdir: full: ok false despite success: " // p1
call remove(p1, ok)
print '(a)', "PASS: mkdir existing"

p1 = "test-filesystem-dir/hello"
call mkdir(p1, ok=ok)
if(.not. is_dir(p1)) error stop "ERROR:test_mkdir: relative: " // p1
if (.not. ok) error stop "ERROR:test_mkdir: relative: ok false despite success: " // p1
call remove(p1, ok)
print '(a)', "PASS: mkdir relative"

end block valgrind

print '(a)', "OK: mkdir"

end program
