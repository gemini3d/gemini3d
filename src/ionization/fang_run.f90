module ionrate
!! these are convenience wrappers for the low-level Fang procedures
!! for if you didn't already have a background atmosphere from the big simulations.

use phys_consts, only: wp
use, intrinsic:: iso_fortran_env, only: sp=>real32
use ionize_fang, only: fang2008, fang2010, gravity_accel, erg2kev
use msis_interface, only : msis_gtd7, msis_gtd8

implicit none (type, external)
private
public :: ionization_fang2008, ionization_fang2010

contains

impure elemental real(wp) function ionization_fang2008(Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, doy, UTsec, &
  msis_version) result(Qtot)

real(wp), intent(in) :: Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, UTsec
integer, intent(in) :: doy, msis_version

real(wp) :: massden_gcm3, meanmass_g
real(wp) :: d(9),T(2), Ap7(7)

Ap7 = Ap

if(msis_version==0) then
  call msis_gtd7(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, d, T, use_meters=.false.)
else
  error stop 'TODO: MSIS 2.x for Fang unit tests'
endif

massden_gcm3 = d(6)  ! [g cm^-3]
meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))

Qtot = fang2008(Q0_erg*erg2kev, E0_keV, T(2), massden_gcm3, meanmass_g, gravity_accel(alt_km))

end function ionization_fang2008


impure elemental real(wp) function ionization_fang2010(Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, doy, UTsec, &
  msis_version) result(Qtot)

real(wp), intent(in) :: Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, UTsec
integer, intent(in) :: doy, msis_version

real(wp) :: massden_gcm3, meanmass_g
real(wp) :: d(9),T(2), Ap7(7)

Ap7 = Ap

if(msis_version==0) then
  call msis_gtd7(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, d, T, use_meters=.false.)
else
  error stop 'TODO: MSIS 2.x for Fang unit tests'
endif

massden_gcm3 = d(6)  ! [g cm^-3]
meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))

Qtot = fang2010(Q0_erg*erg2kev, E0_keV, T(2), massden_gcm3, meanmass_g, gravity_accel(alt_km))

end function ionization_fang2010

end module ionrate
