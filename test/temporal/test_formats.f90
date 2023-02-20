program test_formats
!! This tests various text output formats.
!! It necessarily requires visual inspection, or run from an enclosing program with Regex,
!! such as a Python program configured to test as desired.
use, intrinsic:: iso_fortran_env, only: real32, real64
use phys_consts, only: wp
use timeutils, only: dateinc, utsec2filestem

implicit none (type, external)

character(:), allocatable :: fn
character (25) :: tmp
real(wp) :: UTsec
integer :: ymd(3), i

! print *, 'Single precision lacks adequate precision for dates beyond millisecond'
! print *, utsec2filestem([2015,4,13], 12345.000003_sp)

print *, 'format: easy'
tmp = '20150413_12345.680000'
fn = utsec2filestem([2015,4,13], 12345.678910_wp)
if (fn /= tmp) error stop 'easy: ' // fn

print *, 'format: leading zeros'
tmp = '20150413_00345.680000'
fn = utsec2filestem([2015,4,13],  345.678911_wp)
if (fn /= tmp) error stop 'leading zeros: ' // fn

print *, 'format: leading & trailing zeros'
fn = utsec2filestem([2013,2,20],  60._wp)

if (fn /= '20130220_00060.000000') error stop 'UTsec=60: '//fn

!! Increment time tests
print *, 'format: increment ten_millisecond'
tmp =  '20130220_18000.010000'
ymd = [2013,2,20]
UTsec = 18000._wp
call dateinc(0.01_wp, ymd, UTsec)
fn = utsec2filestem(ymd, UTsec)
if (fn /= tmp) error stop 'format increment ten_millisecond: '//fn

print *, 'format: minute rollover'
tmp =  '20130220_18060.010000'
do i = 1,60
  call dateinc(1.0_wp, ymd, UTsec)
  ! print *, utsec2filestem(ymd, UTsec)
enddo
fn = utsec2filestem(ymd, UTsec)
if (fn /= tmp) error stop 'minute rollover: '//fn

print *, 'format: utsec==86400 corner case'
fn = utsec2filestem([2015,4,13], 86400)

if (fn /= '20150414_00000.000000') error stop 'mismatch 86400 utsec case: '//fn

end program
