program test_copy

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit, int64
use filesystem

implicit none

logical :: ok
integer :: i, u, L1, L2
integer(int64) :: i64, iref
valgrind: block

character(*), parameter :: s1 = "dummy_fortran.txt", s2 = "dummy_fortran.txt.copy"
character(:), allocatable :: t1, t2

t1 = "lots of text today"

open(newunit=u, file=s1, status="replace", action="write", iostat=i)
if(i/=0) error stop "could not open file for writing " // s1
write(u, '(a)') t1
!! format must be "(a)" as "*" makes a blank be prepended to the file
close(u)

!> copy a file
call copy_file(s1, s2, overwrite=.true., ok=ok)
if(.not. is_file(s2)) error stop "did not detect file after copy " // s2
if(.not. ok) error stop "copy_file ok logical is false despite success of copy " // s2

iref = file_size(s1)
i64 = file_size(s2)
if(i64 /= iref) then
  write(stderr, *) "file size mismatch after copy ", s1, s2, iref, i64
  error stop
endif

allocate(character(len(t1)) :: t2)

open(newunit=u, file=s2, status="old", action="read", iostat=i)
if(i/=0) error stop "could not open file for reading " // s2
read(u, '(a)') t2
close(u)

L1 = len_trim(t1)
L2 = len_trim(t2)

if(L1 /= L2 .or. t1 /= t2) then
  write(stderr,*) "file content mismatch after copy: "
  write(stderr,*) s1 // " => "// s2
  write(stderr,*) t1 // " " // t2
  write(stderr,'(a,i0,2x,i0)') "file lengths ", L1, L2
  error stop
endif
!> empty target
! omitted because this fails when ABI shaky e.g. macOS with Clang+Gfortran
! call copy_file(s1, "", status=i)
! if(i==0) error stop "copy_file should fail on empty target"

end block valgrind

print '(a)', "OK: Fortran copyfile"

end program
