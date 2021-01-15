module msis_interface
!! this module allows selecting MSISE00 or MSIS 2.0
!! it is a thin abstraction of the MSIS routines
!! MSISE00 is Fortran 66 style, while MSIS 2.0 is Fortran 90 style
!!
!! We assume MSISE00 is always available, which MSIS 2.0 might not be available.

use, intrinsic :: iso_fortran_env, only : real32, real64
implicit none (type, external)

interface msis_gtd7
  module procedure msis_gtd7_r32, msis_gtd7_r64
end interface msis_gtd7

contains

subroutine msis_gtd7_r32(doy, UTsec, alt_km,  glat, glon, f107a, f107, Ap7, d, T, use_meters, sw25)

external :: meters, gtd7, tselec
integer, intent(in) :: doy
real(real32), intent(in) :: UTsec, alt_km, glat, glon, f107, f107a, Ap7(7)
real(real32), intent(out) :: d(9),T(2)
logical, intent(in) :: use_meters
real(real32), intent(in), optional :: sw25(25)

real(real32) :: stl, sw(25)

stl = UTsec/3600 + glon/15

call meters(use_meters)

sw = 1
if (present(sw25)) sw = sw25
call tselec(sw)

call gtd7(doy, UTsec, alt_km, glat, glon, stl, f107a, f107, Ap7, 48, d, T)

end subroutine msis_gtd7_r32


subroutine msis_gtd7_r64(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, d, T, use_meters, sw25)
!! adds casting to/from real32
external :: meters, gtd7, tselec
integer, intent(in) :: doy
real(real64), intent(in) :: UTsec, alt_km, glat, glon, f107, f107a, Ap7(7)
real(real64), intent(out) :: d(9),T(2)
logical, intent(in) :: use_meters
real(real64), intent(in), optional :: sw25(25)

real(real32) :: sw(25), stl, d32(9), T32(2)

stl = real(UTsec/3600 + glon/15, real32)

call meters(use_meters)

sw = 1
if (present(sw25)) sw = real(sw25, real32)
call tselec(sw)

call gtd7(doy, real(UTsec, real32), real(alt_km, real32), &
  real(glat, real32), real(glon, real32), real(stl, real32), &
  real(f107a, real32), real(f107, real32), real(Ap7, real32), 48, &
  d32, T32)

d = real(d32, real64)
T = real(T32, real64)

end subroutine msis_gtd7_r64

end module
