program pathlib_test

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use pathlib, only : get_filename, mkdir, expanduser, is_absolute, make_absolute, directory_exists

implicit none (type, external)

call test_get_filename()

call test_expanduser_absolute()

call test_directory_exists()


contains


subroutine test_directory_exists()

integer :: i

if(.not.(directory_exists('.'))) error stop "did not detect '.' as directory"

open(newunit=i, file='test-pathlib.h5', status='replace')
close(i)
if((directory_exists('test-pathlib.h5'))) error stop "detected file as directory"
call unlink('test-pathlib.h5')

print *," OK: pathlib: directory_exists"
end subroutine test_directory_exists


subroutine test_get_filename()

character(:), allocatable:: fn
integer :: i
logical :: e

if(get_filename(' ') /= '') error stop 'empty 1'
if(get_filename(' ',' ') /= '') error stop 'empty 2'
!! " " instead of "" to avoid compile-time glitch error with GCC-10 with -Og

call unlink('test-pathlib.h5')

if(len(get_filename('test-pathlib.h5')) > 0) error stop 'not exist full 1'

fn = get_filename('test-pathlib')
if(len(fn) > 0) then
  write(stderr,*) 'ERROR: ',fn, len(fn)
  error stop 'not exist stem 1'
endif

!> touch empty file
open(newunit=i, file='test-pathlib.h5', status='replace')
close(i)
inquire(file='test-pathlib.h5', exist=e)
if(.not.e) error stop 'could not create test-pathlib.h5'

if(get_filename('test-pathlib.h5') /= 'test-pathlib.h5') error stop 'exist full 1'
if(get_filename('test-pathlib') /= 'test-pathlib.h5') error stop 'exist stem 1'

fn = get_filename('.', 'test-pathlib')
if(fn /= './test-pathlib.h5') error stop 'exist stem 2: ' // fn

fn = get_filename('./test-pathlib', 'test-pathlib')
if(fn /= './test-pathlib.h5') error stop 'exist parts 2: ' // fn

fn = get_filename('./test-pathlib.h5', 'test-pathlib')
if(fn /= './test-pathlib.h5') error stop 'exist full 2: ' // fn

call unlink('test-pathlib.h5')

open(newunit=i, file='test-pathlib.nc', status='replace')
close(i)

if(get_filename('test-pathlib.nc') /= 'test-pathlib.nc') error stop 'exist full 1a'
if(get_filename('test-pathlib') /= 'test-pathlib.nc') error stop 'exist stem 1a'
if(get_filename('.', 'test-pathlib') /= './test-pathlib.nc') error stop 'exist stem 2a'
if(get_filename('./test-pathlib', 'test-pathlib') /= './test-pathlib.nc') error stop 'exist parts 2a'

call unlink('test-pathlib.nc')

call mkdir('temp1/temp2')
call unlink('temp1/temp2/test-pathlib.h5')
fn = get_filename('temp1/temp2', 'test-pathlib')
if (fn /= '') error stop 'non-exist dir'
open(newunit=i, file='temp1/temp2/test-pathlib.h5', status='replace')
close(i)
fn = get_filename('temp1/temp2', 'test-pathlib')
if (fn /= 'temp1/temp2/test-pathlib.h5') error stop 'exist dir full 2'

fn = get_filename('./temp1/temp2', 'test-pathlib')
if (fn /= './temp1/temp2/test-pathlib.h5') error stop 'exist dir full 2a'

end subroutine test_get_filename


subroutine test_expanduser_absolute()

character(:), allocatable:: fn
character(16) :: fn2

fn = expanduser("~")
if (fn(1:1) == "/") then
  if (.not.is_absolute("/")) error stop "is_absolute('/') on Unix should be true"
  if (is_absolute("c:/")) error stop "is_absolute('c:/') on Unix should be false"

  fn2 = make_absolute("rel", "/foo")
  if (fn2 /= "/foo/rel") error stop "did not make_absolute Unix /foo/rel, got: " // fn2
else
  if (.not.is_absolute("J:/")) error stop "is_absolute('J:/') on Windows should be true"
  if (.not.is_absolute("j:/")) error stop "is_absolute('j:/') on Windows should be true"
  if (is_absolute("/")) error stop "is_absolute('/') on Windows should be false"

  fn2 = make_absolute("rel", "j:/foo")
  if (fn2 /= "j:/foo/rel") error stop "did not make_absolute Windows j:/foo/rel, got: " // fn2
endif

print *, "OK: pathlib: expanduser,is_absolute"
end subroutine test_expanduser_absolute


subroutine unlink(path)
character(*), intent(in) :: path
integer :: i
logical :: e

inquire(file=path, exist=e)
if (.not.e) return

open(newunit=i, file=path, status='old')
close(i, status='delete')
end subroutine unlink

end program
