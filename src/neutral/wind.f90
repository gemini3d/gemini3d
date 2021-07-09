submodule (neutral) wind
!! https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo/README.txt
use hwm_interface, only : hwm_14, dwm_07
use timeutils, only : ymd2doy

implicit none (type, external)

contains

module procedure neutral_winds

real(wp), dimension(1:size(x%alt,1),1:size(x%alt,2),1:size(x%alt,3)) :: Wmeridional, Wzonal, alt0, v1, v2, v3
integer :: i1,i2,i3, dayOfYear

dayOfYear = ymd2doy(ymd(1), ymd(2), ymd(3))

x1: do i1 = 1,size(x%alt,1)
  x2: do i2 = 1,size(x%alt,2)
    x3: do i3 = 1,size(x%alt,3)
      call hwm_14(dayOfYear, UTsec, &
        alt_km=x%alt(i1,i2,i3), glat=x%glat(i1,i2,i3), glon=x%glon(i1,i2,i3), Ap=Ap, &
        Wmeridional=Wmeridional(i1,i2,i3), Wzonal=Wzonal(i1,i2,i3))
    end do x3
  end do x2
end do x1

alt0 = 0

call rotate_geo2native(vnalt=alt0, vnglat=Wmeridional, vnglon=Wzonal,x=x, vn1=v1, vn2=v2, vn3=v3)

vn1 = vn1 + v1
vn2 = vn2 + v2
vn3 = vn3 + v3

end procedure neutral_winds


end submodule wind
