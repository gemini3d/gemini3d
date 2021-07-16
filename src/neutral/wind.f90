submodule (neutral) wind
!! https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo/README.txt
use hwm_interface, only : hwm_14, dwm_07
use timeutils, only : ymd2doy

implicit none (type, external)

contains

module procedure neutral_winds
  real(wp), dimension(1:size(x%alt,1),1:size(x%alt,2),1:size(x%alt,3)) :: Wmeridional, Wzonal, Walt, v1, v2, v3
  integer :: i1,i2,i3, dayOfYear
  real(wp) :: altnow
  integer :: iinull
  integer :: lx1,lx2,lx3,ix1beg,ix1end
  
  lx1=size(x%alt,1)
  lx2=size(x%alt,2)
  lx3=size(x%alt,3)

  dayOfYear = ymd2doy(ymd(1), ymd(2), ymd(3))
 
  x3: do i3 = 1,lx3
    x2: do i2 = 1,lx2
      x1: do i1 = 1,lx1
        altnow=x%alt(i1,i2,i3)/1.0e3
        if (altnow<0.0) altnow=1.0
        call hwm_14(dayOfYear, UTsec, &
          alt_km=altnow, glat=x%glat(i1,i2,i3), glon=x%glon(i1,i2,i3), Ap=Ap, &
          Wmeridional=Wmeridional(i1,i2,i3), Wzonal=Wzonal(i1,i2,i3))
      end do x1
    end do x2
  end do x3
  
  Walt = 0.0     ! HWM does not provide vertical winds so zero them out
  
  call rotate_geo2native(vnalt=Walt, vnglat=Wmeridional, vnglon=Wzonal,x=x, vn1=v1, vn2=v2, vn3=v3)
  !v1=Walt; v2=Wmeridional; v3=Wzonal;

  !! update module background winds 
  vn1base = v1
  vn2base = v2
  vn3base = v3

  !! zero out background winds at null points
  do iinull=1,x%lnull
    i1=x%inull(iinull,1)
    i2=x%inull(iinull,2)
    i3=x%inull(iinull,3)
    vn1base(i1,i2,i3)=0.0
    vn2base(i1,i2,i3)=0.0
    vn3base(i1,i2,i3)=0.0
  end do

  !! taper the parallel winds toward the null points to avoid artifacts
  do i2=1,lx2
    do i3=1,lx3
      !! identify first and last non-null grid points for this field line
      ix1beg=1
      do while( (.not. x%nullpts(ix1beg,i2,i3)) .and. ix1beg<lx1)     !find the first non-null index for this field line, need to be careful if no null points exist...
        ix1beg=ix1beg+1
      end do
      ix1end=ix1beg
      do while(x%nullpts(ix1end,i2,i3) .and. ix1end<lx1)     !find the first non-null index for this field line
        ix1end=ix1end+1
      end do

      vn1base(ix1beg,i2,i3)=vn1base(ix1beg+1,i2,i3)/2.0
      vn1base(ix1end,i2,i3)=vn1base(ix1end-1,i2,i3)/2.0
    end do
  end do

  !! we really don't resolve mesosphere properly so kill off those winds, these probably don't contribute much to currents???
  !where (x%alt<100e3)
  !  vn1base=0.0
  !  vn2base=0.0
  !  vn3base=0.0
  !end where

  !! force parallel winds to zero to avoid issues...
  !vn1base=0.0     ! it appears to be the case that the parallel drift drives hte model crazy...

  !! update GEMINI wind variables
  call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
end procedure neutral_winds

end submodule wind
