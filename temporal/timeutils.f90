module timeutils

use phys_consts, only: wp, pi

implicit none
private
public :: doy_calc, sza, dateinc

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


pure integer function doy_calc(ymd) result(doy)

integer, intent(in) :: ymd(3)
integer :: i

if ((ymd(3) < 1) .or. (ymd(3) > daysmonth(ymd(1), ymd(2)))) error stop 'impossible day'

doy = 0
do i = 1,ymd(2)-1
  doy = doy + daysmonth(ymd(1), i)
enddo

doy = doy + ymd(3)

end function doy_calc


pure function sza(ymd,UTsec,glat,glon)
!! computes sza in radians
!! CALCULATE SOLAR ZENITH ANGLE OVER A GIVEN GET OF LAT/LON

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: glat,glon

real(wp), dimension(size(glat,1),size(glat,2),size(glat,3)) :: sza

real(wp) :: doy,soldecrad
real(wp), dimension(size(glat,1),size(glat,2),size(glat,3)) :: lonrad,LThrs,latrad,hrang

!> SOLAR DECLINATION ANGLE
doy=doy_calc(ymd)
soldecrad=-23.44_wp*cos(2._wp*pi/365._wp*(doy+10._wp))*pi/180._wp;

!> HOUR ANGLE
lonrad=glon*pi/180._wp
where (lonrad>pi)
  lonrad=lonrad-2._wp*pi
end where
where (lonrad<-2._wp*pi)
  lonrad=lonrad+2._wp*pi
end where
LThrs=UTsec/3600._wp+lonrad/(pi/12._wp)
hrang=(12._wp-LThrs)*(pi/12._wp)

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
if (utsec < 0.) error stop 'negative UTsec, simulation should go forward in time only!'
if (dtsec < 0.) error stop 'negative dtsec, simulation should go forward in time only!'
if (dtsec > 86400.) error stop 'excessively large dtsec, simulation step should be small enough!'

UTsec = UTsec + dtsec
if (UTsec >= 86400._wp) then
  UTsec = modulo(UTsec, 86400._wp)
  day = day+1          !roll the day over

!> month rollover
  if (day > daysmonth(year, month)) then
    day=1
    monthinc=1
  else
    monthinc=0
  end if

  month = month + monthinc

  if (month>12) then    !< roll the year over
    month=1
    year=year+1
  end if
end if

ymd(1)=year; ymd(2)=month; ymd(3)=day;    !replace input array with new date

end subroutine dateinc

end module timeutils
