submodule (neutral) wind
!! https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo/README.txt
use hwm_interface, only : hwm_14, dwm_07
use timeutils, only : ymd2doy

implicit none (type, external)

contains
!>  This procedure makes the call to HWM for horizontal winds.  It loads those winds into the neutral module background
!     arrays but does not directly assign them to any "output" variables as it used to - this avoids the need, strictly
!     spreaking to have this procedure call another procedure that is dependend on mpi. At the same time it does require
!     the "main" program to make an additional call to assign the background (and any perturbations) to variables used in that
!     program for neutral parameters.
module procedure neutral_winds
  real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3) :: Wmeridional, Wzonal, Walt
  integer :: i1,i2,i3, dayOfYear
  real(wp) :: altnow,glonnow,glatnow
  integer :: iinull
  integer :: lx1,lx2,lx3

  lx1=x%lx1
  lx2=x%lx2
  lx3=x%lx3

  dayOfYear = ymd2doy(ymd(1), ymd(2), ymd(3))

  x3: do i3 = 1,lx3
    x2: do i2 = 1,lx2
      x1: do i1 = 1,lx1
        glonnow=x%glon(i1,i2,i3)
        glatnow=x%glat(i1,i2,i3)
        altnow=x%alt(i1,i2,i3)/1.0e3
        if (altnow<0.0) altnow=1.0
        call hwm_14(dayOfYear, UTsec, &
          alt_km=altnow, glat=glatnow, glon=glonnow, Ap=Ap, &
          Wmeridional=Wmeridional(i1,i2,i3), Wzonal=Wzonal(i1,i2,i3))
      end do x1
    end do x2
  end do x3

  Walt = 0.0     ! HWM does not provide vertical winds so zero them out

  !print*, 'Rotating atmospheric information from HWM14'
  call rotate_geo2native(vnalt=Walt, vnglat=Wmeridional, vnglon=Wzonal,x=x, atmos=atmos, flagBG=.true.)
  !v1=Walt; v2=Wmeridional; v3=Wzonal;

  !! zero out background winds at null points
  do iinull=1,x%lnull
    i1=x%inull(iinull,1)
    i2=x%inull(iinull,2)
    i3=x%inull(iinull,3)
    atmos%vn1base(i1,i2,i3)=0.0
    atmos%vn2base(i1,i2,i3)=0.0
    atmos%vn3base(i1,i2,i3)=0.0
  end do

  !! taper winds according to altitude.  If this is not done there seems to be an issue where poorly resolved
  !!  drifts in the lower E-region cause stability problems.  Generally speaking, it's not too bad to omit field-
  !!  aligned drifts in the E-region since most of the dynamical behavior there is driven by the field-perp winds
  !!  (which are retained).  That being said, this could have implications, e.g. for spE modeling so perhaps should
  !!  be revisited in the future.
  do i2=1,lx2
    do i3=1,lx3
      atmos%vn1base(1:lx1,i2,i3)=atmos%vn1base(1:lx1,i2,i3)*(0.5 + 0.5*tanh((x%alt(1:lx1,i2,i3)-150e3)/10e3))
    end do
  end do

  !! update GEMINI wind variables
  !call neutral_wind_update(vn1,vn2,vn3,v2grid,v3grid)
end procedure neutral_winds

end submodule wind
