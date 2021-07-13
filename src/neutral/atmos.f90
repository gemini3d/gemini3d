submodule (neutral) atmos

use, intrinsic :: iso_fortran_env, only: sp => real32, stderr=>error_unit

use timeutils, only : ymd2doy
use msis_interface, only : msis_gtd7, msis_gtd8

implicit none (type, external)

contains

!! subroutine neutral_atmos(ymd,UTsecd,glat,glon,alt,activ,nn,Tn)
!! CALL NRL-MSISE-00 AND ORGANIZE THE RESULTS.  APPEND
!! OTHER AUXILIARY NEUTRAL DENSITY DATA USED BY MAIN CODE
module procedure neutral_atmos
  integer :: ix1,ix2,ix3,lx1,lx2,lx3
  integer :: doy
  real(wp) :: ap(7),ap3
  real(wp) :: altnow
  real(wp) :: d(9),t(2)
    !   real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: nnow
  !    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: altalt    !an alternate altitude variable which fixes below ground values to 1km
  
  lx1=size(alt,1)
  lx2=size(alt,2)
  lx3=size(alt,3)
  
  !! CONVERT DATE INFO INTO EXPECTED FORM AND KIND
  ap = activ(3)
  ap3 = activ(3)
  doy = ymd2doy(year=ymd(1), month=ymd(2), day=ymd(3))
  ap(2)=ap3   !superfluous for now
  
  !> ITERATED LAT, LON, ALT DATA
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        altnow= alt(ix1,ix2,ix3)/1000
        if (altnow < 0) then
          altnow = 1     !so that MSIS does not get called with below ground values and so that we set them to something sensible that won't mess up the conductance calculations
        end if
  
        if(msis_version == 0) then
          call msis_gtd7(doy=doy, UTsec=UTsecd, &
            alt_km=altnow, glat=glat(ix1,ix2,ix3), glon=glon(ix1,ix2,ix3), &
            f107a=activ(1), f107=activ(2), ap7=ap, &
            d=d, T=t, use_meters=.true.)
        elseif(msis_version == 20) then
          call msis_gtd8(doy=doy, UTsec=UTsecd, &
            alt_km=altnow, glat=glat(ix1,ix2,ix3), glon=glon(ix1,ix2,ix3), &
            f107a=activ(1), f107=activ(2), ap7=ap, &
            Dn=d, Tn=t)
        else
          write(stderr,*) 'ERROR:neutral_atmos: unknown msis version',msis_version,' expected 0 (MSISE00) or 20 (MSIS 2.0)'
          error stop
        end if
  
        nnmsis(ix1,ix2,ix3,1)= d(2)
        nnmsis(ix1,ix2,ix3,2)= d(3)
        nnmsis(ix1,ix2,ix3,3)= d(4)
        nnmsis(ix1,ix2,ix3,4)= d(7)
        nnmsis(ix1,ix2,ix3,5)= d(8)
  
        Tnmsis(ix1,ix2,ix3)= t(2)
        nnmsis(ix1,ix2,ix3,6)=0.4_wp*exp(-3700/Tnmsis(ix1,ix2,ix3))*nnmsis(ix1,ix2,ix3,3)+ &
                            5e-7_wp*nnmsis(ix1,ix2,ix3,1)   !Mitra, 1968
      end do
    end do
  end do
  
  !> if HWM selected call it similar to MSIS above, if not zero all out (default)
  !vn1base = 0; vn2base = 0; vn3base = 0
  !! see winds submodule for HWM handling
  
  !> Update current state with new background and existing perturbations, if used
  call neutral_denstemp_update(nn,Tn)
end procedure neutral_atmos

end submodule atmos
