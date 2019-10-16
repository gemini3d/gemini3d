!! test fang2008 implementation
!! trying to reproduce data of Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
use, intrinsic:: iso_fortran_env, only: dp=>real64, sp=>real32
use ionize_fang, only: fang2008
implicit none

integer :: i, iyd
real(sp), parameter :: alt_km(77) = [(real(i, sp), i=20,400,5)]
real(dp), parameter :: E0_keV(5) = [0.1_dp, 1._dp, 10._dp, 100._dp, 1000._dp]
real(dp), parameter :: Q0_erg = 1._dp
real(dp), parameter :: erg2kev = 624150648._dp
real(dp) :: massden_gcm3, meanmass_g, g_ms2
real(sp) :: f107, f107a, ap(7), glat, glon, sec, stl, d(9),T(2), sw(25)
real(dp), parameter :: Re_m = 6371.0e3_dp, Me_kg = 5.9722e24_dp, Gconst=6.67408e-11_dp
real(dp) :: Qtot(size(alt_km), size(E0_keV))

!! MSIS setup
iyd = 00001
sec = 0
glat = 60
glon = 0
f107 = 50
f107a = f107
ap = 5
STL=SEC/3600+GLON/15
call meters(.false.)

sw = 1
call tselec(sw)
call tretrv(sw)
print '(A,/,25F3.0)', 'MSIS SW:', SW

do i = 1, size(alt_km)
  call gtd7(iyd,sec,alt_km(i),glat,glon,stl,f107a,f107,ap,48,d,T)
  massden_gcm3 = d(6)  ! [g cm^-3]
  meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))
  g_ms2 = -1 * Gconst * Me_kg / (Re_m + alt_km(i)*1000)**2  ! [m s^-2]

  ! print *, massden_gcm3, meanmass_g
  ! massden_gcm3 ~ 1e-5..1e-15
  ! meanmass_g ~ 3e-26

  Qtot(i, :) = fang2008(Q0_erg*erg2kev, E0_keV, real(T(2), dp), massden_gcm3, meanmass_g, g_ms2)
enddo

print '(A)', 'alt[km]/E0[keV]    0.1            1.0             10.           100.            1000.'
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot(i,:)
enddo


end program