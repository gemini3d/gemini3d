module timeutils
use phys_consts, only: wp
use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
implicit none
private
public :: doy_calc, sza, dateinc

real(wp), parameter :: pi = 4._wp*atan(1._wp)

contains

elemental integer function daysmonth(year, month) result(days)

integer, intent(in) :: year, month
integer :: monthdays(12)

if ((year < 1600) .or. (year > 2500)) error stop 'is year specified correctly?'
if ((month < 1) .or. (month > 12)) error stop 'impossible month'

monthdays = [31,28,31,20,31,30,31,31,30,31,30,31]

if (mod(year, 4)==0 .and. mod(year, 100)/=0 .or. mod(year, 400) == 0) monthdays(2)=29

days = monthdays(month)

end function daysmonth


elemental integer function doy_calc(year, month, day) result(doy)

integer, intent(in) :: year, month, day
integer :: i

if ((day < 1) .or. (day > daysmonth(year, month))) error stop 'impossible day'

doy = 0
do i = 1, month-1
  doy = doy + daysmonth(year, i)
enddo

doy = doy + day

end function doy_calc


elemental function sza(year, month, day, UTsec,glat,glon)
!! computes sza in radians
!! CALCULATE SOLAR ZENITH ANGLE OVER A GIVEN GET OF LAT/LON

integer, intent(in) :: year, month, day
real(wp), intent(in) :: UTsec
real(wp), intent(in) :: glat,glon
real(wp) :: sza

real(wp), parameter :: tau = 2._wp*pi

real(wp) :: doy,soldecrad
real(wp) :: lonrad,LThrs,latrad,hrang

!> SOLAR DECLINATION ANGLE
doy=doy_calc(year, month, day)
soldecrad=-23.44_wp*cos(tau/365._wp*(doy+10)) * pi/180

!> HOUR ANGLE
lonrad=glon*pi/180
lonrad = modulo(lonrad, 2*pi)

LThrs=UTsec/3600._wp+lonrad/(pi/12._wp)
hrang=(12-LThrs)*(pi/12._wp)

!> SOLAR ZENITH ANGLE
latrad=glat*pi/180._wp
sza=acos(sin(soldecrad)*sin(latrad)+cos(soldecrad)*cos(latrad)*cos(hrang))

end function sza


subroutine dateinc(dtsec, ymd, UTsec)

real(wp), intent(in) :: dtsec
!! seconds to increment
integer, intent(inout) :: ymd(3)
!! year, month, day of month
real(wp), intent(inout) :: UTsec
!! seconds since midnight UTC

integer :: year,month,day
integer :: monthinc          !< we incremented the month

year=ymd(1); month=ymd(2); day=ymd(3);

if ((day < 1) .or. (day > 31)) error stop 'impossible day'
if (utsec < 0) error stop 'negative UTsec, simulation should go forward in time only!'
if (dtsec < 0) error stop 'negative dtsec, simulation should go forward in time only!'
if (dtsec > 86400) error stop 'excessively large dtsec, simulation step should be small enough!'

UTsec = UTsec + dtsec
if (UTsec >= 86400) then
  UTsec = modulo(UTsec, 86400._wp)
  day = day+1          !roll the day over

!> month rollover
  if (day > daysmonth(year, month)) then
    day=1
    month = month + 1
  end if

  if (month>12) then    !< roll the year over
    month=1
    year=year+1
  end if
end if

ymd(1)=year; ymd(2)=month; ymd(3)=day;    !replace input array with new date

end subroutine dateinc

end module timeutils
