module timeutils
use phys_consts, only: wp
use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64, int32, int64
implicit none
private
public :: doy_calc, sza, dateinc, utsec2filestem, date_filename, day_wrap

real(wp), parameter :: pi = 4._wp*atan(1._wp)

contains

elemental integer function daysmonth(year, month) result(days)

integer, intent(in) :: year, month
integer :: monthdays(12)

if ((year < 1600) .or. (year > 2500)) error stop 'is year specified correctly?'
if ((month < 1) .or. (month > 12)) error stop 'impossible month'

monthdays = [31,28,31,30,31,30,31,31,30,31,30,31]

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


pure subroutine dateinc(dtsec, ymd, UTsec)
!! increment datetime by dtsec

real(wp), intent(in) :: dtsec
!! seconds to increment
integer, intent(inout) :: ymd(3)
!! year, month, day of month
real(wp), intent(inout) :: UTsec
!! seconds since midnight UTC

integer :: year,month,day
integer :: monthinc          !< we incremented the month

year=ymd(1); month=ymd(2); day=ymd(3);

if (day < 1) error stop 'temporal:timeutils:dateinc(): days are positive integers'
if (utsec < 0) error stop 'negative UTsec, simulation should go forward in time only!'
if (dtsec < 0) error stop 'negative dtsec, simulation should go forward in time only!'
if (dtsec > 86400) error stop 'excessively large dtsec > 86400, simulation step should be small enough!'

UTsec = UTsec + dtsec
do while (UTsec >= 86400)
  UTsec = UTsec - 86400._wp
  day = day+1
  call day_wrap(year, month, day)
end do

ymd(1)=year; ymd(2)=month; ymd(3)=day;    !replace input array with new date

end subroutine dateinc


recursive pure subroutine day_wrap(year, month, day)
!! increment date if needed, according to day
!! that is, if day is beyond month, increment month and year if needed
integer, intent(inout) :: year, month, day

if (month < 1 .or. day < 1) error stop 'day_wrap: months and days must be positive'

!> wrap months
do while (month > 12)
  month = month - 12
  year = year + 1
end do

!> wrap days
do while (day > daysmonth(year, month))
  day = day - daysmonth(year, month)
  month = month + 1
  call day_wrap(year, month, day)
end do

end subroutine day_wrap


pure function date_filename(outdir,ymd,UTsec)  result(filename)
!! GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
character(:), allocatable :: filename


!> assemble
filename = outdir // '/' // utsec2filestem(ymd, UTsec)

end function date_filename


pure character(21) function utsec2filestem(ymd, UTsec) result(fn)
!! file stem is exactly 21 characters long, per Matt Z's de facto spec.
!! FIXME: until we go to integer UTsec (microsec) we round to nearest microsecond
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)
integer(int64) :: usec, seconds, frac
character(12) :: sec_str
integer(int64), parameter :: million = 1000000_int64
integer :: year, month, day
real(dp) :: UTsectmp

year = ymd(1)
month = ymd(2)
day = ymd(3)

select type(UTsec)
  type is (real(sp))
    usec = nint(UTsec * million, int64)
  type is (real(dp))
    usec = nint(UTsec * million, int64)
  type is (integer(int32))
    usec = int(UTsec, int64) * million
  type is (integer(int64))
    usec = UTsec * million
  class default
    error stop "io/formats.f90:utsec2filestem unknown UTsec type"
end select

seconds = usec / million
if (seconds < 0 .or. seconds > 86400) error stop 'io/formats.f90:utsec2filestem did NOT satisfy 0 <= seconds < 86400'
if (seconds == 86400) then
  !> FIXME: This corner case is from not using integers for microseconds
  ! write(stderr,*) 'utsec2filestem: FIXME: patching UTsec=86400 to next day midnight'
  day = day+1
  seconds = 0
  usec = 0
  call day_wrap(year, month, day)
endif
frac = usec - seconds * million

write(sec_str, '(I5.5, A1, I6.6)') seconds, '.', frac

write(fn,'(i4,I2.2,I2.2,a13)') year, month, day, '_' // sec_str

end function utsec2filestem

end module timeutils
