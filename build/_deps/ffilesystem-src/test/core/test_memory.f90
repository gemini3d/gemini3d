program test_owner

use, intrinsic :: iso_fortran_env

use filesystem

implicit none

integer(int64) :: f, t
logical :: ok = .true.

f = free_memory()
t = total_sys_memory()

if (f <= 0) then
  write(error_unit, "(a, i0)") 'free_memory() returned non-positive value', f
  ok = .false.
end if
if (t <= 0) then
  write(error_unit, "(a, i0)") 'total_sys_memory() returned non-positive value', t
  ok = .false.
end if
if (f > t) then
  write(error_unit, "(a, i0, a, i0)") 'free_memory() returned value larger than total_sys_memory: ', f, ' > ', t
  ok = .false.
end if

if (.not. ok) error stop 'memory test failed.'

print '(a,i0)', 'Free memory: ', f
print '(a,i0)', 'Total system memory: ', t

print '(a)', 'OK: memory test passed.'

end program
