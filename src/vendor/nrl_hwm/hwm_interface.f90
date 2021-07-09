module hwm_interface
!! abstract interface to NRL HWM horizontal wind model

use, intrinsic :: iso_fortran_env, only : real32, real64

implicit none (type, external)

external :: hwm14, dwm07


interface hwm_14
  procedure :: hwm_14_r64, hwm_14_r32
end interface

interface dwm_07
  procedure :: dwm_07_r64, dwm_07_r32
end interface

private
public :: hwm_14, dwm_07

contains

subroutine hwm_14_r64(dayOfYear, UTsec, alt_km, glat, glon, Ap, Wmeridional, Wzonal)
!! Parameters
!! ----------
!!
!! dayOfYear
!! UTsec
!! alt_km
!! glat : geodetic latitude(deg)
!! glon : geodetic longitude(deg)
!! ap : current 3hr ap index
!!
!! Returns
!! -------
!! w(1) = meridional wind (m/sec + northward)
!! w(2) = zonal wind (m/sec + eastward)
!!
!! Like MSIS, HWM does not use the year in iyd. Just give the day of year.

integer, intent(in) :: dayOfYear
real(real64), intent(in) :: UTsec, alt_km, glat, glon, Ap
real(real64), intent(out) :: Wmeridional, Wzonal

real(real32) :: Ap2(2), W(2), dummy

Ap2(2) = real(Ap, real32)

call hwm14(dayOfYear, real(UTsec, real32), &
  real(alt_km, real32), real(glat, real32), real(glon, real32), &
  dummy, dummy, dummy, Ap2, W)

Wmeridional = real(W(1), real64)
Wzonal = real(W(2), real64)

end subroutine hwm_14_r64


subroutine hwm_14_r32(dayOfYear, UTsec, alt_km, glat, glon, Ap, Wmeridional, Wzonal)
!! Parameters
!! ----------
!!
!! dayOfYear
!! UTsec
!! alt_km
!! glat : geodetic latitude(deg)
!! glon : geodetic longitude(deg)
!! ap : current 3hr ap index
!!
!! Returns
!! -------
!! w(1) = meridional wind (m/sec + northward)
!! w(2) = zonal wind (m/sec + eastward)
!!
!! Like MSIS, HWM does not use the year in iyd. Just give the day of year.

integer, intent(in) :: dayOfYear
real(real32), intent(in) :: UTsec, alt_km, glat, glon, Ap
real(real32), intent(out) :: Wmeridional, Wzonal

real(real32) :: Ap2(2), dummy, W2(2)

Ap2(2) = Ap

call hwm14(dayOfYear, UTsec, alt_km, glat, glon, dummy, dummy, dummy, Ap2, W2)

Wmeridional = W2(1)
Wzonal = W2(2)

end subroutine hwm_14_r32


subroutine dwm_07_r64(dayOfYear, UTsec, alt_km, glat, glon, Ap, DW2)

integer, intent(in) :: dayOfYear
real(real64), intent(in) ::  UTsec, alt_km, glat, glon, Ap
real(real64), intent(out) :: DW2(2)

real(real32) :: DW(2), Ap2(2)

Ap2(2) = real(Ap, real32)

call dwm07(dayOfYear,real(UTsec, real32), &
  real(alt_km, real32), real(glat, real32), real(glon, real32), &
  Ap2, DW)

DW2 = real(DW, real64)

end subroutine dwm_07_r64

subroutine dwm_07_r32(dayOfYear, UTsec, alt_km, glat, glon, Ap, DW2)

integer, intent(in) :: dayOfYear
real(real32), intent(in) ::  UTsec, alt_km, glat, glon, Ap
real(real32), intent(out) :: DW2(2)

real(real32) :: Ap2(2)

Ap2(2) = Ap

call dwm07(dayOfYear,UTSEC,ALT_km,GLAT,GLON,AP2,DW2)

end subroutine dwm_07_r32


end module hwm_interface
