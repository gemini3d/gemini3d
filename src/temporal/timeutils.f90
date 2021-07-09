module timeutils
use phys_consts, only: wp
use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64, int32, int64

implicit none (type, external)
private
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


pure function date_filename(outdir, ymd, UTsec) result(filename)
!! GENERATE A FILENAME stem OUT OF A GIVEN DATE/TIME
!! (does not include suffix like .h5)

character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
class(*), intent(in) :: UTsec

character(:), allocatable :: filename
character(len=21) :: stem

stem = utsec2filestem(ymd, UTsec)

filename = outdir // '/' // stem

end function date_filename


pure character(21) function utsec2filestem(ymd, UTsec) result(fn)
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

frac = millisec*1000 - seconds * 1000000  !< microseconds

write(sec_str, '(I5.5, A1, I6.6)') seconds, '.', frac

write(fn,'(i4,I2.2,I2.2,a13)') year, month, day, '_' // sec_str

end function utsec2filestem


pure subroutine find_lastdate(ymd0,UTsec0,ymdtarget,UTsectarget,cadence,ymd,UTsec)

!> Compute the last date before the target date based on a start date and date cadence.  The
!  file is assumed to exist and programmer needs to check for existence outside this
!  procedure.
! FIXME:  this will occasionally hang and may need to be re-examined at some point...

integer, dimension(3), intent(in) :: ymd0
real(wp), intent(in) :: UTsec0
integer, dimension(3), intent(in) :: ymdtarget
real(wp), intent(in) :: UTsectarget
real(wp), intent(in) :: cadence
integer, dimension(3), intent(out) :: ymd
real(wp), intent(out) :: UTsec

integer, dimension(3) :: ymdnext
real(wp) :: UTsecnext
logical :: flagend


ymd=ymd0
UTsec=UTsec0
ymdnext=ymd0
UTsecnext=UTsec0
flagend=ymdnext(1)>=ymdtarget(1) .and. ymdnext(2)>=ymdtarget(2) .and. ymdnext(3)>=ymdtarget(3) .and. UTsecnext>UTsectarget &
          .or. ymdnext(1)>ymdtarget(1) .or. ymdnext(2)>ymdtarget(2) .or. ymdnext(3)>ymdtarget(3) ! in case the first time step is the last before target
do while ( .not.(flagend) )
  ymd=ymdnext
  UTsec=UTsecnext
  call dateinc(cadence,ymdnext,UTsecnext)
  flagend=ymdnext(1)>=ymdtarget(1) .and. ymdnext(2)>=ymdtarget(2) .and. ymdnext(3)>=ymdtarget(3) .and. UTsecnext>UTsectarget &
            .or. ymdnext(1)>ymdtarget(1) .or. ymdnext(2)>ymdtarget(2) .or. ymdnext(3)>ymdtarget(3)
end do
! When the loops exits ymd,UTsec have the date of the last output file before the given target time OR the first output file in the case that the target date is before the begin date...

end subroutine find_lastdate


pure function find_time_elapsed(ymdstart,UTsecstart,ymdend,UTsecend,dt) result(telapsed)

! finds the amount of time that has elapsed between a given start and end date, using a given dt increment
! the resulting elapsed time will be the smallest multiple of dt that is >= true elapsed time

!! inputs
integer, dimension(3), intent(in) :: ymdstart,ymdend
real(wp), intent(in) :: UTsecstart,UTsecend
reaL(wp), intent(in) :: dt

integer, dimension(3) :: ymdnow
real(wp) :: UTsecnow

real(wp) :: telapsed  !! output

telapsed=0._wp; ymdnow=ymdstart; UTsecnow=UTsecstart;
do while (.not. (all(ymdend==ymdnow) .and. UTsecnow>UTsecend) )
  call dateinc(dt,ymdnow,UTsecnow)
  telapsed=telapsed+dt
end do
telapsed=telapsed-dt

end function find_time_elapsed


end module timeutils
