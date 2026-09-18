program test_canon

use, intrinsic:: iso_fortran_env, only : stderr=>error_unit

use filesystem

implicit none

valgrind : block

character(:), allocatable :: p1, p2

integer :: L1

! -- current directory  -- old MacOS doesn't handle "." or ".." alone
p1 = realpath(".")
L1 = len_trim(p1)
if(L1 == 0) error stop "FAILED: realpath '.'"

p2 = realpath(get_cwd())
if(p1 /= p2) then
  write(stderr, '(a)')  "FAILED: realpath '.' failed: " // p1 // " /= " // p2
  error stop
endif

print *, "OK: current dir = ", p1


! -- relative dir
p1 = realpath("..")

p2 = parent(p2)
if (p1 /= p2) then
  write(stderr,'(a)') 'ERROR:realpath:relative: up dir not realpath ' // p1 // " /= " // p2
  error stop
end if
print *, 'OK: realpath(..) = ', p1

end block valgrind


end program
