program test

use, intrinsic:: iso_fortran_env, only : stderr=>error_unit

use filesystem

implicit none


valgrind : block

logical :: ok

character(:), allocatable :: old_cwd, cwd

character(:), allocatable :: s1, s2

s1 = getarg(0)
print '(a)', 'program name: '// s1

s2 = parent(parent(s1))
print '(a)', 'parent(parent(' // s1 // ')) = '// s2

old_cwd = get_cwd()

print '(a)', 'current working directory: '// old_cwd

ok = set_cwd(s2)

if (.not. ok) error stop "chdir failed"

if(set_cwd("")) error stop "chdir should fail on empty string"

cwd = get_cwd()
print '(a)', 'New working directory: '// cwd

if (.not. same_file(cwd, s2)) error stop "chdir failed: " // s2 // " /= " // cwd

end block valgrind

print '(a)', "OK: chdir"

end program
