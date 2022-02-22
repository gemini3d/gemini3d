module hwm_interface
!! when HWM14 is disabled, we simply output zeros.

use, intrinsic :: iso_fortran_env, only : real32, real64

implicit none (type, external)

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

integer, intent(in) :: dayOfYear
real(real64), intent(in) :: UTsec, alt_km, glat, glon, Ap
real(real64), intent(out) :: Wmeridional, Wzonal

Wmeridional = 0
Wzonal = 0

end subroutine hwm_14_r64


subroutine hwm_14_r32(dayOfYear, UTsec, alt_km, glat, glon, Ap, Wmeridional, Wzonal)

integer, intent(in) :: dayOfYear
real(real32), intent(in) :: UTsec, alt_km, glat, glon, Ap
real(real32), intent(out) :: Wmeridional, Wzonal

Wmeridional = 0
Wzonal = 0

end subroutine hwm_14_r32


subroutine dwm_07_r64(dayOfYear, UTsec, alt_km, glat, glon, Ap, DW2)

integer, intent(in) :: dayOfYear
real(real64), intent(in) ::  UTsec, alt_km, glat, glon, Ap
real(real64), intent(out) :: DW2(2)

DW2 = 0

end subroutine dwm_07_r64


subroutine dwm_07_r32(dayOfYear, UTsec, alt_km, glat, glon, Ap, DW2)

integer, intent(in) :: dayOfYear
real(real32), intent(in) ::  UTsec, alt_km, glat, glon, Ap
real(real32), intent(out) :: DW2(2)

DW2 = 0

end subroutine dwm_07_r32


end module hwm_interface
