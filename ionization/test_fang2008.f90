!! test fang2008 implementation
!! trying to reproduce data of Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
!!
!! This could readily be turned into a subroutine. For now we use it as a registration case.
use, intrinsic:: iso_fortran_env, only: dp=>real64, sp=>real32
use ionize_fang, only: fang2008, gravity_accel, erg2kev
implicit none

integer :: i, iyd
real(sp), parameter :: alt_km(*) = [(real(i, sp), i=20,400,5)]
real(dp), parameter :: E0_keV(*) = [0.1_dp, 1._dp, 10._dp, 100._dp, 1000._dp], &
                       Q0_erg = 1._dp
real(dp) :: massden_gcm3, meanmass_g
real(sp) :: f107, f107a, ap(7), glat, glon, sec, stl, d(9),T(2), sw(25)
real(dp), allocatable :: Qtot(:,:)

allocate(Qtot(size(alt_km), size(E0_keV)))

!! MSIS setup
iyd = 00001
sec = 0
glat = 60
glon = 0
f107 = 50
f107a = f107
ap(:) = 5
STL=SEC/3600+GLON/15
call meters(.false.)

sw(:) = 1
call tselec(sw)
call tretrv(sw)
print '(A,/,25F3.0)', 'MSIS SW:', SW

do i = 1, size(alt_km)
  call gtd7(iyd,sec,alt_km(i),glat,glon,stl,f107a,f107,ap,48,d,T)
  massden_gcm3 = d(6)  ! [g cm^-3]
  meanmass_g = massden_gcm3 / (sum(d(1:5)) + sum(d(7:8)))

  ! print *, massden_gcm3, meanmass_g
  ! massden_gcm3 ~ 1e-5..1e-15
  ! meanmass_g ~ 3e-26

  Qtot(i, :) = fang2008(Q0_erg*erg2kev, E0_keV, real(T(2), dp), massden_gcm3, meanmass_g, real(gravity_accel(alt_km(i)), dp))
enddo

print '(A)', 'alt[km]/E0[keV]    0.1            1.0             10.           100.            1000.'
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot(i,:)
enddo


end program