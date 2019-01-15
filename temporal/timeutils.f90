module timeutils

use phys_consts, only: wp, pi

implicit none
private
public :: doy_calc, sza, dateinc

contains

pure integer function doy_calc(ymd)

integer, intent(in) :: ymd(3)

integer :: monthday(12)

monthday=[31,28,31,20,31,30,31,31,30,31,30,31]
if (mod(ymd(1),4)==0) then
  monthday(2)=29  !< leap year adjustment
end if
if(ymd(2)==1) then
  doy_calc=ymd(3)
else
  doy_calc=sum(monthday(1:ymd(2)-1))+ymd(3) 
end if

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


subroutine dateinc(dt,ymd,UTsec)
!! TODO: the leap year computation is incomplete: 100,400 year rules

real(wp), intent(in) :: dt
integer, intent(inout) :: ymd(3)
real(wp), intent(inout) :: UTsec

integer :: year,month,day
integer :: monthinc          !let's us know whether we incremented the month
integer :: daymax     !number of days February should have for given year

year=ymd(1); month=ymd(2); day=ymd(3);


UTsec=UTsec+dt
if (UTsec>=86400._wp) then
  UTsec=modulo(UTsec,86400._wp)
  day=day+1          !roll the day over

  select case(month)     !check whether the month needs to be rolled over
    case(4,6,9,11)    !months having 30 days
      daymax=30
    case(2)    !Feb.
      if (mod(year,4)==0) then     !leap year
        daymax=29
      else          !normal year
        daymax=28
      end if 
    case default    !eveyone else has 31 days
      daymax=30
  end select

  if (day>daymax) then    !roll the month over
    day=1
    month=month+1
    monthinc=1
  else
    monthinc=0
  end if

  if (monthinc==1 .and. month>12) then    !roll the year over
    month=1
    year=year+1
  end if
end if

ymd(1)=year; ymd(2)=month; ymd(3)=day;    !replace input array with new date

end subroutine dateinc

end module timeutils
