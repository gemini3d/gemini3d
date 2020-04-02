module ionrate
!! these are convenience wrappers for the low-level Fang procedures
!! for if you didn't already have a background atmosphere from the big simulations.

use phys_consts, only: wp
use, intrinsic:: iso_fortran_env, only: sp=>real32
use ionize_fang, only: fang2008, fang2010, gravity_accel, erg2kev

implicit none (external)
private
public :: ionization_fang2008, ionization_fang2010

contains

impure elemental real(wp) function ionization_fang2008(Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, doy, UTsec) result(Qtot)

real(wp), intent(in) :: Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, UTsec
integer, intent(in) :: doy

real(wp) :: massden_gcm3, meanmass_g
real(sp) :: stl, d(9),T(2), sw(25), Ap7(7)

external :: meters, gtd7, tselec

!! MSIS setup
Ap7 = real(Ap, sp)
stl = real(UTsec/3600 + glon/15, sp)
call meters(.false.)
sw = 1
call tselec(sw)
call gtd7(doy, real(UTsec, sp), real(alt_km, sp), real(glat, sp), real(glon, sp), stl, &
          real(f107a, sp), real(f107, sp), Ap7, 48, d, T)
massden_gcm3 = d(6)  ! [g cm^-3]
meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))

Qtot = fang2008(Q0_erg*erg2kev, E0_keV, real(T(2), wp), massden_gcm3, meanmass_g, real(gravity_accel(alt_km), wp))

end function ionization_fang2008


impure elemental real(wp) function ionization_fang2010(Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, doy, UTsec) result(Qtot)

real(wp), intent(in) :: Q0_erg, E0_keV, alt_km, f107, f107a, Ap, glat, glon, UTsec
integer, intent(in) :: doy

real(wp) :: massden_gcm3, meanmass_g
real(sp) :: stl, d(9),T(2), sw(25), Ap7(7)

external :: meters, gtd7, tselec

!! MSIS setup
Ap7 = real(Ap, sp)
stl = real(UTsec/3600 + glon/15, sp)
call meters(.false.)
sw = 1
call tselec(sw)
call gtd7(doy, real(UTsec, sp), real(alt_km, sp), real(glat, sp), real(glon, sp), stl, &
          real(f107a, sp), real(f107, sp), Ap7, 48, d, T)
massden_gcm3 = d(6)  ! [g cm^-3]
meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))

Qtot = fang2010(Q0_erg*erg2kev, E0_keV, real(T(2), wp), massden_gcm3, meanmass_g, real(gravity_accel(alt_km), wp))

end function ionization_fang2010

end module