submodule (neutral) atmos

use, intrinsic :: iso_fortran_env, only: sp => real32, stderr=>error_unit

use timeutils, only : ymd2doy
use msis_interface, only : msis_gtd7, msis_gtd8

implicit none (type, external)

contains
!>  This procedure makes the call to MSIS for density/temperature.  It loads those data into the neutral module background
!     arrays but does not directly assign them to any "output" variables as it used to - this avoids the need, strictly
!     spreaking to have this procedure call another procedure that is dependend on mpi.  At the same time it does require
!     the "main" program to make an additional call to assign the background (and any perturbations) to variables used in that
!     program for neutral parameters.
module procedure neutral_atmos
  integer :: ix1,ix2,ix3,lx1,lx2,lx3
  integer :: doy
  real(wp) :: ap(7),ap3
  real(wp) :: altnow,glonnow,glatnow
  real(wp) :: d(9),t(2)
    !   real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: nnow
  !    real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: altalt    !an alternate altitude variable which fixes below ground values to 1km

  lx1=size(atmos%nnmsis,1)
  lx2=size(atmos%nnmsis,2)
  lx3=size(atmos%nnmsis,3)

  !! CONVERT DATE INFO INTO EXPECTED FORM AND KIND
  ap = activ(3)
  ap3 = activ(3)
  doy = ymd2doy(year=ymd(1), month=ymd(2), day=ymd(3))
  ap(2)=ap3   !superfluous for now

  !> ITERATED LAT, LON, ALT DATA, note that we keep calling MSIS even for periodic grids but since this isn't triggered every
  !    time step we avoid generating more voluminous (but efficient) code
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        glonnow=glon(ix1,ix2,ix3)
        glatnow=glat(ix1,ix2,ix3)
        altnow= alt(ix1,ix2,ix3)/1000

        if (altnow < 0.0) then
          altnow = 1.0     !so that MSIS does not get called with below ground values and so that we set them to something sensible that won't mess up the conductance calculations
        end if

        if(msis_version == 0) then
          !! MSISE00
          call msis_gtd7(doy=doy, UTsec=UTsecd, &
            alt_km=altnow, glat=glatnow, glon=glonnow, &
            f107a=activ(1), f107=activ(2), ap7=ap, &
            d=d, T=t, use_meters=.true.)
        else
          !! MSIS 2.x
          call msis_gtd8(doy=doy, UTsec=UTsecd, &
            alt_km=altnow, glat=glatnow, glon=glonnow, &
            f107a=activ(1), f107=activ(2), ap7=ap, &
            Dn=d, Tn=t)
        end if

        atmos%nnmsis(ix1,ix2,ix3,1)= d(2)
        atmos%nnmsis(ix1,ix2,ix3,2)= d(3)
        atmos%nnmsis(ix1,ix2,ix3,3)= d(4)
        atmos%nnmsis(ix1,ix2,ix3,4)= d(7)
        atmos%nnmsis(ix1,ix2,ix3,5)= d(8)

        atmos%Tnmsis(ix1,ix2,ix3)= t(2)
        atmos%nnmsis(ix1,ix2,ix3,6)=0.4_wp*exp(-3700/atmos%Tnmsis(ix1,ix2,ix3))*atmos%nnmsis(ix1,ix2,ix3,3)+ &
                            5e-7_wp*atmos%nnmsis(ix1,ix2,ix3,1)   !Mitra, 1968
      end do
    end do
  end do

  !> if HWM selected call it similar to MSIS above, if not zero all out (default)
  !vn1base = 0; vn2base = 0; vn3base = 0
  !! see winds submodule for HWM handling

  !> Update current state with new background and existing perturbations, if used
  !call neutral_denstemp_update(nn,Tn)
end procedure neutral_atmos

end submodule atmos
