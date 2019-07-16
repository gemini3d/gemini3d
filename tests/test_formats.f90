!! This tests various text output formats.
!! It necessarily requires visual inspection, or run from an enclosing program with Regex,
!! such as a Python program configured to test as desired.
use, intrinsic:: iso_fortran_env, only: sp=>real32, dp=>real64, stderr=>error_unit
use date_formats, only: utsec2filename
use timeutils, only: dateinc

implicit none

character(:), allocatable :: fn
real(dp) :: UTsec, dtinc
integer :: ymd(3), i

! print *, 'Single precision lacks adequate precision for dates beyond millisecond'
! print *, utsec2filename([2015,4,13], 12345.000003_sp)

fn = utsec2filename([2015,4,13], 12345.678910_dp)

if (fn /= '20150413_12345.678910.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch real64'
endif


fn = utsec2filename([2015,4,13],  345.678911_dp)

if (fn /= '20150413_00345.678911.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch leading zeros'
endif


fn = utsec2filename([2013,2,20],  60._dp)

if (fn /= '20130220_00060.000000.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch UTsec=60'
endif

!! Increment time tests
ymd = [2013,2,20]
UTsec = 18000._dp
call dateinc(1e-6_dp, ymd, UTsec)
fn = utsec2filename(ymd, UTsec)
if (fn /= '20130220_18000.000001.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch output filename'
endif

do i = 1,60
  call dateinc(1._dp, ymd, UTsec)
  ! print *, utsec2filename(ymd, UTsec)
enddo
fn = utsec2filename(ymd, UTsec)
if (fn /= '20130220_18060.000001.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch leading zeros'
endif


end program
