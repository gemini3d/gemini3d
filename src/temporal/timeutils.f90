module timeutils
use phys_consts, only: wp
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64, int32, int64

implicit none (type, external)
private
public :: elapsed_seconds, shift_datetime
public :: ymd2doy, sza, dateinc, utsec2filestem, date_filename, day_wrap, find_lastdate, find_time_elapsed

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


elemental integer function ymd2doy(year, month, day) result(doy)

integer, intent(in) :: year, month, day
integer :: i

if ((day < 1) .or. (day > daysmonth(year, month))) error stop 'impossible day'

doy = 0
do i = 1, month-1
  doy = doy + daysmonth(year, i)
enddo

doy = doy + day

end function ymd2doy


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
doy = ymd2doy(year, month, day)
soldecrad = -23.44_wp*cos(tau/365._wp*(doy+10)) * pi/180

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

year=ymd(1); month=ymd(2); day=ymd(3);

if (ymd2doy(year,month,day)<1) error stop 'timeutils: invalid date'
if (.not.all(ieee_is_finite([dtsec,UTsec]))) error stop 'timeutils: nonfinite time'
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


pure function date_filename(outdir, ymd, UTsec)
!! GENERATE A FILENAME stem OUT OF A GIVEN DATE/TIME
!! (does not include suffix like .h5)

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
character(:), allocatable :: date_filename

character(len=21) :: stem

stem = utsec2filestem(ymd, UTsec)

date_filename = outdir // '/' // stem

end function date_filename


pure character(21) function utsec2filestem(ymd, UTsec)
!! file stem is exactly 21 characters long, per Matt Z's de facto spec.
!! we keep microsecond dummy precision filenames to be legacy compatible
!! true filename resolution is 10 milliseconds due to real32 7 digits of precision vis 86400 second day.
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec
!! UTC second: real [0.0 .. 86400.0)

character(12) :: sec_str
integer :: year, month, day, seconds, millisec, frac

year = ymd(1)
month = ymd(2)
day = ymd(3)

select type(UTsec)
  type is (real(sp))
    !! round to nearest ten milliseconds
    millisec = nint(UTsec * 100) * 10
  type is (real(dp))
    !! round to nearest ten milliseconds
    millisec = nint(UTsec * 100) * 10
  type is (integer(int32))
    millisec = UTsec * 1000
  type is (integer(int64))
    millisec = int(UTsec) * 1000
  class default
    error stop "timeutils.f90:utsec2filestem unknown UTsec type"
end select

seconds = millisec / 1000 !< truncate fractional second
if (seconds < 0 .or. seconds > 86400) error stop 'timeutils.f90::utsec2filestem did NOT satisfy 0 <= seconds < 86400'
if (seconds == 86400) then
  !> FIXME: This corner case is from not using integers for microseconds
  ! write(stderr,*) 'utsec2filestem: FIXME: patching UTsec=86400 to next day midnight'
  day = day+1
  seconds = 0
  millisec = 0
  call day_wrap(year, month, day)
endif

frac = modulo(millisec,1000) * 1000  ! microseconds without signed integer overflow

write(sec_str, '(I5.5, A1, I6.6)') seconds, '.', frac

write(utsec2filestem, '(i4,I2.2,I2.2,a13)') year, month, day, '_' // sec_str

end function utsec2filestem


!> Gregorian day number, with input validation through ymd2doy.
pure integer(int64) function day_number(ymd) result(n)
integer, intent(in) :: ymd(3)
integer(int64) :: y
y=int(ymd(1)-1,int64)
n=365_int64*y+y/4-y/100+y/400+ymd2doy(ymd(1),ymd(2),ymd(3))
end function day_number

!> Signed elapsed seconds; no iteration over time steps and no loss of date at midnight.
pure real(wp) function elapsed_seconds(ymd0,ut0,ymd1,ut1) result(dt)
integer, intent(in) :: ymd0(3),ymd1(3)
real(wp), intent(in) :: ut0,ut1
if (.not.all(ieee_is_finite([ut0,ut1]))) error stop 'timeutils: nonfinite UTC'
if (min(ut0,ut1)<0 .or. max(ut0,ut1)>=86400) error stop 'timeutils: UTC outside [0,86400)'
dt=real(day_number(ymd1)-day_number(ymd0),wp)*86400._wp+ut1-ut0
end function elapsed_seconds

!> Shift a normalized UTC date in either direction, including month/year boundaries.
pure subroutine shift_datetime(offset,ymd,utsec)
real(wp), intent(in) :: offset
integer, intent(inout) :: ymd(3)
real(wp), intent(inout) :: utsec
integer :: ndays,i
real(wp) :: total
if (.not.all(ieee_is_finite([offset,utsec]))) error stop 'timeutils: nonfinite offset'
if (abs(offset)>366._wp*86400*900) error stop 'timeutils: offset exceeds supported calendar'
if (ymd2doy(ymd(1),ymd(2),ymd(3))<1) error stop 'timeutils: invalid date'
if (utsec<0 .or. utsec>=86400) error stop 'timeutils: UTC outside [0,86400)'
total=utsec+offset
ndays=floor(total/86400._wp)
utsec=total-real(ndays,wp)*86400._wp
do i=1,abs(ndays)
  if (ndays>0) then
    ymd(3)=ymd(3)+1
    call day_wrap(ymd(1),ymd(2),ymd(3))
  else
    ymd(3)=ymd(3)-1
    if (ymd(3)==0) then
      ymd(2)=ymd(2)-1
      if (ymd(2)==0) then
        ymd(1)=ymd(1)-1; ymd(2)=12
      endif
      ymd(3)=daysmonth(ymd(1),ymd(2))
    endif
  endif
enddo
if (ymd2doy(ymd(1),ymd(2),ymd(3))<1) error stop 'timeutils: date outside supported range'
end subroutine shift_datetime

!> Latest cadence timestamp at or before target; return start when target precedes it.
pure subroutine find_lastdate(ymd0,UTsec0,ymdtarget,UTsectarget,cadence,ymd,UTsec)
integer, intent(in) :: ymd0(3),ymdtarget(3)
real(wp), intent(in) :: UTsec0,UTsectarget,cadence
integer, intent(out) :: ymd(3)
real(wp), intent(out) :: UTsec
real(wp) :: dt,offset
if (.not.ieee_is_finite(cadence)) error stop 'timeutils: nonfinite cadence'
if (cadence<=0 .or. cadence>86400) error stop 'timeutils: cadence must be in (0,86400]'
dt=max(0._wp,elapsed_seconds(ymd0,UTsec0,ymdtarget,UTsectarget))
offset=real(floor(dt/cadence,kind=int64),wp)*cadence
ymd=ymd0; UTsec=UTsec0
call shift_datetime(offset,ymd,UTsec)
end subroutine find_lastdate

!> Greatest nonnegative multiple of dt not exceeding elapsed time (legacy semantics).
pure real(wp) function find_time_elapsed(ymdstart,UTsecstart,ymdend,UTsecend,dt) result(telapsed)
integer, intent(in) :: ymdstart(3),ymdend(3)
real(wp), intent(in) :: UTsecstart,UTsecend,dt
real(wp) :: actual
if (.not.ieee_is_finite(dt)) error stop 'timeutils: nonfinite cadence'
if (dt<=0 .or. dt>86400) error stop 'timeutils: cadence must be in (0,86400]'
actual=elapsed_seconds(ymdstart,UTsecstart,ymdend,UTsecend)
if (actual<0) error stop 'timeutils: end precedes start'
telapsed=real(floor(actual/dt,kind=int64),wp)*dt
end function find_time_elapsed

end module timeutils
