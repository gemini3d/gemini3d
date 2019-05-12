!! test day rollover
use phys_consts, only: wp
use timeutils, only: dateinc, doy_calc

implicit none

integer :: ymd(3)
real(wp) :: dtsec, UTsec

!> leap year tests
if (doy_calc(2015,1,1) /= 1) error stop 'doy_calc day 1'

if (doy_calc(2015,2,28) /= 59) error stop 'doy_calc day 59'

if (doy_calc(2012,2,29) /= 60) error stop 'doy_calc leap day 60'

if (doy_calc(2000,2,29) /= 60) error stop 'doy_calc millenium leap day 60'

if (doy_calc(1900,3,1) /= 60) error stop 'doy_calc 1900 NO leap day 60'

if (doy_calc(2100,3,1) /= 60) error stop 'doy_calc 2100 NO leap day 60'


!> date rollover tests
ymd = [1900,2,28]
dtsec = 0.
utsec = 86400.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [1900,3,1])) error stop 'intrinsic rollover'

ymd = [1900,2,28]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [1900,3,1])) error stop '1900 NO leap year rollover'

ymd = [2012,2,28]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2012,2,29])) error stop '2012 leap year rollover'


ymd = [1999,12,31]
dtsec = 1.
utsec = 86399.
call dateinc(dtsec, ymd, utsec)
if (.not.all(ymd == [2000,1,1]) .and. utsec<=1e-3) error stop '2000 leap year rollover'

end program
