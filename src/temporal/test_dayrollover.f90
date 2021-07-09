program test_dayrollover
!! test day rollover
use, intrinsic :: iso_fortran_env, only: real32, real64, stderr=>error_unit
use phys_consts, only: wp
use timeutils, only: dateinc, ymd2doy, day_wrap

implicit none (type, external)

integer :: ymd(3), year, month, day
real(wp) :: dtsec, UTsec

!> leap year tests
if (ymd2doy(2015,1,1) /= 1) error stop 'ymd2doy day 1'

if (ymd2doy(2015,2,28) /= 59) error stop 'ymd2doy day 59'

if (ymd2doy(2012,2,29) /= 60) error stop 'ymd2doy leap day 60'

if (ymd2doy(2000,2,29) /= 60) error stop 'ymd2doy millennium leap day 60'

if (ymd2doy(1900,3,1) /= 60) error stop 'ymd2doy 1900 NO leap day 60'

if (ymd2doy(2100,3,1) /= 60) error stop 'ymd2doy 2100 NO leap day 60'

!> day rollover checks
year = 2000
month = 2
day = 30
call day_wrap(year, month, day)
if (year /= 2000 .or. month /= 3 .or. day /= 1) error stop 'day_wrap: single day leap fail'

year = 2000
month = 2
day = 31
call day_wrap(year, month, day)
if (year /= 2000 .or. month /= 3 .or. day /= 2) error stop 'day_wrap: two day leap fail'

year = 1999
month = 14
day = 31
call day_wrap(year, month, day)
if (year /= 2000 .or. month /= 3 .or. day /= 2) error stop 'day_wrap: 14 month + two day leap fail'


!> date rollover tests
ymd = [1900,2,28]
dtsec = 0.
utsec = 86400.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [1900,3,1]) .and. abs(utsec) < epsilon(1.)) error stop 'intrinsic rollover'

ymd = [1900,2,28]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [1900,3,1]) .and. abs(utsec) < epsilon(1.)) error stop '1900 NO leap year rollover'

ymd = [2012,2,28]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2012,2,29]) .and. abs(utsec) < epsilon(1.)) error stop '2012 leap year rollover'


print *,'rollover 1999-12-31T24:00:00'
ymd = [1999,12,31]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2000,1,1]) .and. abs(utsec) < epsilon(1.)) error stop '2000 century rollover'

ymd = [2000,2,28]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2000,2,29]) .and. abs(utsec) < epsilon(1.)) error stop '2000 leap year rollover'

!> multi-day rollover
ymd = [2000,2,28]
dtsec = 1._wp
utsec = 172799._wp
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2000,3,1]) .and. abs(utsec) < epsilon(1.)) error stop '2000 leap year rollover'

ymd = [2000,2,1]
dtsec = 1.
utsec = 2592000._wp
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2000,3,1]) .and. abs(utsec) < epsilon(1.)) error stop '2000 leap year + month rollover'

if (wp == real64) then
  print *, 'rollover: 366 days leap year'
  ymd = [2000,2,1]
  dtsec = 1.
  utsec = 31622400._wp
  call dateinc(dtsec, ymd, utsec)
  if (.not.all(ymd == [2001,3,1]) .and. abs(utsec) < epsilon(1.)) then
    write(stderr, *) 'FAILED: expected 2001-03-01 but got', ymd
    error stop 'FAILED: 2000 leap year 366 day rollover'
  endif
else
  write(stderr, *) 'SKIPPED: 366 day rollover test: real32 not being precise enough'
endif

print *, 'rollover: 10 years + leap year'
ymd = [2000,2,1]
dtsec = 1.
utsec = 316224000._wp
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2010,2,8]) .and. abs(utsec) < epsilon(1.)) error stop '10 years day rollover'


print *, 'rollover: 100 years + leap year'
ymd = [2000,2,1]
dtsec = 1.
utsec = 3162240000._wp
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2100,4,17]) .and. abs(utsec) < epsilon(1.)) error stop '100 years day rollover'


end program
