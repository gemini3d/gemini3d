program pathlib_test

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use pathlib, only : get_filename, mkdir

implicit none (type, external)

character(:), allocatable :: fn
integer :: i
logical :: e


!> test get_filename

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
if(fn /= './test-pathlib.h5') then
  write(stderr,*) fn
  error stop 'exist stem 2'
endif

fn = get_filename('./test-pathlib', 'test-pathlib')
if(fn /= './test-pathlib.h5') then
  write(stderr,*) fn
  error stop 'exist parts 2'
endif

fn = get_filename('./test-pathlib.h5', 'test-pathlib')
if(fn /= './test-pathlib.h5') then
  write(stderr,*) fn
  error stop 'exist full 2'
endif

call unlink('test-pathlib.h5')

open(newunit=i, file='test-pathlib.nc', status='replace')
close(i)

if(get_filename('test-pathlib.nc') /= 'test-pathlib.nc') error stop 'exist full 1a'
if(get_filename('test-pathlib') /= 'test-pathlib.nc') error stop 'exist stem 1a'
if(get_filename('.', 'test-pathlib') /= './test-pathlib.nc') error stop 'exist stem 2a'
if(get_filename('./test-pathlib', 'test-pathlib') /= './test-pathlib.nc') error stop 'exist parts 2a'

call unlink('test-pathlib.nc')

i = mkdir('temp1/temp2')
call unlink('temp1/temp2/test-pathlib.h5')
fn = get_filename('temp1/temp2', 'test-pathlib')
if (fn /= '') error stop 'non-exist dir'
open(newunit=i, file='temp1/temp2/test-pathlib.h5', status='replace')
close(i)
fn = get_filename('temp1/temp2', 'test-pathlib')
if (fn /= 'temp1/temp2/test-pathlib.h5') error stop 'exist dir full 2'

fn = get_filename('./temp1/temp2', 'test-pathlib')
if (fn /= './temp1/temp2/test-pathlib.h5') error stop 'exist dir full 2a'

print *, 'OK: pathlib'

contains

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
