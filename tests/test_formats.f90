!! This tests various text output formats.
!! It necessarily requires visual inspection, or run from an enclosing program with Regex,
!! such as a Python program configured to test as desired.
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp
use timeutils, only: dateinc, utsec2filename

implicit none

character(:), allocatable :: fn
real(wp) :: UTsec, dtinc
integer :: ymd(3), i

! print *, 'Single precision lacks adequate precision for dates beyond millisecond'
! print *, utsec2filename([2015,4,13], 12345.000003_sp)

print *, 'format: easy'
fn = utsec2filename([2015,4,13], 12345.678910_wp)

if (fn /= '20150413_12345.678910.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch real64'
endif


print *, 'format: leading zeros'
fn = utsec2filename([2015,4,13],  345.678911_wp)

if (fn /= '20150413_00345.678911.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch leading zeros'
endif


print *, 'format: leading & trailing zeros'
fn = utsec2filename([2013,2,20],  60._wp)

if (fn /= '20130220_00060.000000.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch UTsec=60'
endif

!! Increment time tests
print *, 'format: increment microsecond'
ymd = [2013,2,20]
UTsec = 18000._wp
call dateinc(1e-6_wp, ymd, UTsec)
fn = utsec2filename(ymd, UTsec)
if (fn /= '20130220_18000.000001.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch output filename'
endif

print *, 'format: minute rollover'
do i = 1,60
  call dateinc(1.0_wp, ymd, UTsec)
  ! print *, utsec2filename(ymd, UTsec)
enddo
fn = utsec2filename(ymd, UTsec)
if (fn /= '20130220_18060.000001.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch leading zeros'
endif


print *, 'format: utsec==86400 corner case'
fn = utsec2filename([2015,4,13], 86400)

if (fn /= '20150414_00000.000000.dat') then
  write(stderr,*) 'wrong output: '//fn
  error stop 'mismatch 86400 utsec case'
endif

end program
