submodule (neutral) atmos

use, intrinsic :: iso_fortran_env, only: sp => real32

use timeutils, only : doy_calc

implicit none (type, external)

contains

module procedure neutral_atmos
!! subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,nn,Tn)
!! CALL NRL-MSISE-00 AND ORGANIZE THE RESULTS.  APPEND
!! OTHER AUXILIARY NEUTRAL DENSITY DATA USED BY MAIN CODE

integer :: ix1,ix2,ix3,lx1,lx2,lx3

integer :: iyd,mass=48
integer :: dom,month,year,doy,yearshort
real(sp) :: sec,f107a,f107,ap(7),stl,ap3
real(sp) :: altnow,latnow,lonnow
real(sp) :: d(9),t(2)

external :: meters, gtd7

!   real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: nnow
!    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: altalt    !an alternate altitude variable which fixes below ground values to 1km


lx1=size(alt,1)
lx2=size(alt,2)
lx3=size(alt,3)


!! CONVERT DATE INFO INTO EXPECTED FORM AND KIND
f107a=real(activ(1),sp)
f107=real(activ(2),sp)
ap=real(activ(3),sp)
ap3=real(activ(3),sp)
dom=ymd(3)
month=ymd(2)
year=ymd(1)
doy=doy_calc(year, month, dom)
yearshort=mod(year,100)
iyd=yearshort*1000+doy
sec=floor(UTsecd)

ap(2)=ap3   !superfluous for now


!> ITERATED LAT, LON, ALT DATA
call meters(.true.)    !switch to mksa units

do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      altnow=real(alt(ix1,ix2,ix3)/1d3,sp)
      if (altnow<0.0) then
        altnow = 1._sp     !so that MSIS does not get called with below ground values and so that we set them to something sensible that won't mess up the conductance calculations
      end if

      latnow=real(glat(ix1,ix2,ix3),sp)
      lonnow=real(glon(ix1,ix2,ix3),sp)

      stl=sec/3600.0+lonnow/15.0
      call gtd7(iyd,sec,altnow,latnow,lonnow,stl,f107a,f107,ap,mass,d,t)

      nnmsis(ix1,ix2,ix3,1)=real(d(2),wp)
      nnmsis(ix1,ix2,ix3,2)=real(d(3),wp)
      nnmsis(ix1,ix2,ix3,3)=real(d(4),wp)
      nnmsis(ix1,ix2,ix3,4)=real(d(7),wp)
      nnmsis(ix1,ix2,ix3,5)=real(d(8),wp)

      Tnmsis(ix1,ix2,ix3)=real(t(2),wp)
      nnmsis(ix1,ix2,ix3,6)=4d-1*exp(-3700.0/Tnmsis(ix1,ix2,ix3))*nnmsis(ix1,ix2,ix3,3)+ &
                          5d-7*nnmsis(ix1,ix2,ix3,1)   !Mitra, 1968
    end do
  end do
end do


!> Could be replaced with HWM at some point, but for now make these zero
vn1base = 0d0; vn2base = 0d0; vn3base = 0d0


!> Update current state with new background and existing perturbations
call neutral_update(nn,Tn,vn1,vn2,vn3)

!> subtract off grid drift speed (needs to be set to zero if not lagrangian grid)
vn2=vn2-v2grid
vn3=vn3-v3grid

end procedure neutral_atmos

end submodule atmos
