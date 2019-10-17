!! Reproduces data of:
!! * Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
!! * Figure 2 in Fang 2010

use, intrinsic:: iso_fortran_env, only: sp=>real32
use phys_consts, only: wp
use ionrate, only: ionization_fang2008, ionization_fang2010
implicit none

integer :: i, iyd
real(wp), parameter :: alt_km(*) = [(real(i, wp), i=20,400,2)]
real(wp), parameter :: E0_keV(*) = [0.1, 1., 10., 100., 1000.], &
                       Q0_erg = 1.
real(wp), parameter :: f107=50, f107a=50, Ap=5, glat=60, glon=0, utsec=0
integer, parameter :: doy=1
real(wp) :: massden_gcm3, meanmass_g
real(wp), allocatable :: Qtot08(:,:), Qtot10(:,:)

allocate(Qtot08(size(alt_km), size(E0_keV)), Qtot10(size(alt_km), size(E0_keV)))


do i = 1, size(alt_km)
  Qtot08(i, :) = ionization_fang2008(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec)
  Qtot10(i, :) = ionization_fang2010(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec)
enddo

print '(A)', 'alt[km]/E0[keV]    0.1            1.0             10.           100.            1000.'
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot08(i,:)
enddo

print '(/,A)', 'alt[km]/Emono[keV]  0.1            1.0             10.           100.            1000.'
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot10(i,:)
enddo


end program