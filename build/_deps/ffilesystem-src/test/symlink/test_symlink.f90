program test_symlink

use, intrinsic:: iso_fortran_env, only : stderr=>error_unit
use filesystem

implicit none

valgrind: block

logical :: ok

character(:), allocatable :: tgt, rtgt, link, linko, tgt_dir, link_dir

tgt_dir = parent(getarg(0))

!> create and verify target file
tgt = realpath(tgt_dir // filesep() // "test_symlink_fortran.txt")

print '(a)', "target file: " // tgt

if (is_file(tgt)) then
  print *, "deleting old target " // tgt
  call remove(tgt)
end if

call touch(tgt, ok)
if (.not. ok) then
  write(stderr, "(a)") "ERROR: touch("//tgt//") failed"
  error stop 77
end if
print '(a)', "created target file " // tgt

if(.not. is_file(tgt)) then
  write(stderr, "(a)") "ERROR: is_file("//tgt//") should be true for existing regular file target"
  error stop
end if

link = join(tgt_dir, "test_symlink_fortran.link")
linko = join(tgt_dir, "test_oo_fortran.link")
link_dir = join(tgt_dir, "my_link_fortran.dir")

! print *, "TRACE:test_symlink: target: " // tgt
! print *, "TRACE:test_symlink: link: " // link

call create_symlink(tgt, "", ok)
if (ok) error stop "ERROR: create_symlink() should fail with empty link"
print '(a)', "PASSED: create_symlink: empty link"
if(is_symlink(tgt)) then
  write(stderr, '(a)') "is_symlink() should be false for non-symlink file: " // tgt
  error stop
end if

call create_symlink("", link, ok)
if (ok) error stop "ERROR: create_symlink() should fail with empty target"
print '(a)', "PASSED: create_symlink: empty target"

if (is_symlink(link)) then
  print '(a)', "deleting old symlink " // link
  call remove(link)
end if
call create_symlink(tgt, link)
print '(a)', "PASSED: create_symlink " // link

!> read_symlink
rtgt = read_symlink(link)
if(len_trim(rtgt) /= len_trim(tgt)) then
  write(stderr, '(a)') "read_symlink() failed: " // rtgt // " /= " // tgt
  error stop
end if
print '(a)', "PASSED: read_symlink " // rtgt // " == " // tgt

!> read_symlink non-symlink
rtgt = read_symlink(tgt)
if (len_trim(rtgt) > 0) then
  write(stderr, '(a)') "read_symlink() should return empty string for non-symlink file: " // rtgt
  error stop
end if

!> read_symlink non-existent
rtgt = read_symlink("not-exist-file")
if (len_trim(rtgt) > 0) then
  write(stderr, '(a)') "read_symlink() should return empty string for non-existent file: " // rtgt
  error stop
end if



if (is_symlink(linko)) then
  print *, "deleting old symlink " // linko
  call remove(linko)
end if
call create_symlink(tgt, linko)
print '(a)', "PASSED: created symlink " // linko

!> directory symlinks
if (is_symlink(link_dir)) then
  print *, "deleting old symlink " // link_dir
  call remove(link_dir)
end if
call create_symlink(tgt_dir, link_dir)

!> checks
! call create_symlink("", "")  !< this error stops

!> file symlinks
if(is_symlink(tgt)) then
  write(stderr, '(a)') "ERROR: is_symlink("//tgt//") should be false for non-symlink target"
  error stop
end if

if(.not. is_symlink(link)) then
  write(stderr, '(a)') "ERROR: is_symlink("//link//") should be true"
  error stop
end if

if(.not. is_file(link)) then
  write(stderr, "(a)") "ERROR: is_file("//link//") should be true for existing regular file target " // tgt
  error stop
end if

print '(a)', "PASSED: test_symlink: file"

!> directory symlinks
if(is_symlink(tgt_dir)) error stop "is_symlink() should be false for non-symlink dir"
if(.not. is_dir(link_dir)) then
  write(stderr, '(a)') "ERROR: is_dir("//link_dir//") should be true for existing regular dir"
  error stop
end if

if(.not. is_symlink(link_dir)) then
  write(stderr, '(a)') "ERROR: is_symlink() should be true for symlink dir: " // link_dir
  error stop
end if

end block valgrind

print *, "OK: filesystem symbolic links"

end program
