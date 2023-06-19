module ionize_fang

use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
use phys_consts, only: wp, kb

implicit none (type, external)
private
public :: fang2008, fang2010, fang2010_spectrum, gravity_accel, erg2kev

real(wp), parameter :: deps = 0.035_wp
real(wp), parameter :: erg2kev = 624150648._wp
!! keV, kinetic energy lost per ion-electron pair produced

contains


elemental real(wp) function fang2010_spectrum(Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g_ms2, diff_num_flux, kappa &
, bimax_frac, E0_char_keV) result(qtot)
!! composite spectrum of isotropically precipitating monoenergetic (100 eV to 1 MeV) electrons
!! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GL045406

real(wp), intent(in) :: Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g_ms2, kappa, bimax_frac, E0_char_keV
integer, intent(in) :: diff_num_flux

real(wp) :: y, H_cm, f, dQ0_keV, phi_keV, Ebin_keV, dEbin_keV, bin_lb, bin_ub
integer :: i, j, k
real(wp), dimension(8) :: C
character(12) :: E0_str
integer, parameter :: lbins = 64

real(wp), parameter :: P(8,4) = reshape( &
[1.24616_wp,     1.45903_wp,    -2.42269e-1_wp,  5.95459e-2_wp, &
 2.23976_wp,    -4.22918e-7_wp,  1.36458e-2_wp,  2.53332e-3_wp, &
 1.41754_wp,     1.44597e-1_wp,  1.70433e-2_wp,  6.39717e-4_wp, &
 2.48775e-1_wp, -1.50890e-1_wp,  6.30894e-9_wp,  1.23707e-3_wp, &
-4.65119e-1_wp, -1.05081e-1_wp, -8.95701e-2_wp,  1.22450e-2_wp, &
 3.86019e-1_wp,  1.75430e-3_wp, -7.42960e-4_wp,  4.60881e-4_wp, &
-6.45454e-1_wp,  8.49555e-4_wp, -4.28581e-2_wp, -2.99302e-3_wp, &
 9.48930e-1_wp,  1.97385e-1_wp, -2.50660e-3_wp, -2.06938e-3_wp], shape(P), order=[2,1])

 real(wp), parameter :: phi_Evans_keV(64) = &
 [1.99356e8_wp,	1.82982e8_wp,	1.70856e8_wp,	1.60908e8_wp, &
  1.52842e8_wp,	1.46431e8_wp,	1.41496e8_wp,	1.39092e8_wp, &
  1.36728e8_wp,	1.36728e8_wp,	1.36728e8_wp,	1.34404e8_wp, &
  1.34404e8_wp,	1.34404e8_wp,	1.34404e8_wp,	1.34404e8_wp, &
  1.34404e8_wp,	1.34404e8_wp,	1.34404e8_wp,	1.34404e8_wp, &
  1.34404e8_wp,	1.34404e8_wp,	1.34404e8_wp,	1.36728e8_wp, &
  1.39092e8_wp,	1.42714e8_wp,	1.46431e8_wp,	1.47692e8_wp, &
  1.51538e8_wp,	1.54158e8_wp,	1.59534e8_wp,	1.63689e8_wp, &
  1.67953e8_wp,	1.73810e8_wp,	1.78337e8_wp,	1.86145e8_wp, &
  1.17517e9_wp,	9.98575e8_wp,	8.26979e8_wp,	6.73232e8_wp, &
  5.34156e8_wp,	4.09528e8_wp,	3.22155e8_wp,	2.36630e8_wp, &
  1.65098e8_wp,	1.17182e8_wp,	0.76341e8_wp,	0.50594e8_wp, &
  0.31849e8_wp,	0.15336e8_wp,	0.08840e8_wp,	0.04054e8_wp, &
  0.02064e8_wp,	0.00832e8_wp,	0.00374e8_wp,	0.00132e8_wp, &
  4.66259e4_wp,	1.29068e4_wp,	3.54600e3_wp,	8.60021e2_wp, &
  1.84025e2_wp,	3.47254e1_wp,	5.77661e0_wp,	7.51004e-1_wp]
 
if (E0_keV < 0.1_wp .or. E0_keV > 1000) then
  write(E0_str,'(F12.4)') E0_keV
  error stop 'ionize_fang:fang2010_spectrum: valid E0 range from 100 eV .. 1 MeV: E0 (keV) : ' // E0_str
endif

!! scale height
!! Equation (2)
H_cm = 100 * kb * Tn / (meanmass_g/1000) / abs(g_ms2)

qtot = 0
do k=1,lbins
  !! energy bin limits
  bin_lb = -1 ! lower limit from Fang 2008, 2010
  select case (diff_num_flux)
  case (0) ! Maxwellian
    bin_ub = log10(20*E0_keV) ! gives 1.12e-7 x max(phi)
  case (3) ! accelerated Maxwellian
    bin_ub = log10(20*E0_keV)
  case (4) ! Evans, D. S. (1974) 2 keV acc., 0.8 keV temp
    bin_ub = log10(10*2.0_wp)
  case default
    bin_ub = 3 ! upper limit from Fang 2008, 2010
  end select

  !! log bins, midpoint rule
  Ebin_keV =  (10**(bin_lb+(bin_ub-bin_lb)*(k)/(lbins-1)) + 10**(bin_lb+(bin_ub-bin_lb)*(k-1)/(lbins-1)))/2 
  dEbin_keV = (10**(bin_lb+(bin_ub-bin_lb)*(k)/(lbins-1)) - 10**(bin_lb+(bin_ub-bin_lb)*(k-1)/(lbins-1)))

  !! normalized atmospheric column mass
  !! Equation (1)
  y = 2/Ebin_keV * (massden_gcm3 * H_cm / 6e-6_wp)**0.7_wp

  !! Equation (5)
  C = 0
  do i=1,size(P,1)
    do j=1,size(P,2)
      C(i) = C(i) + P(i,j) * log(Ebin_keV)**(j-1)
    end do
  end do
  C = exp(C)

  !! Equation (4)
  !! Energy deposition "f"
  f = C(1)*y**C(2)*exp(-1*C(3)*y**C(4)) + C(5)*y**C(6) * exp(-1*C(7)*y**C(8))

  !! user choice of differential energy flux
  select case (diff_num_flux)
  case (0) ! Maxwellian
    phi_keV = (Q0_keV/(2*E0_keV**3)) * Ebin_keV * exp(-1*Ebin_keV/E0_keV) ! 1/keV/s/cm^2
  case (1) ! kappa distribution
    if (kappa <= 2) then
      error stop 'ERROR:ionize_fang:fang2010_spectrum: for finite <E>, kappa must be greater than 2'
    end if
    phi_keV = (Q0_keV/(2*E0_keV**3)) * (gamma(kappa+1)/(gamma(kappa-1/2)*(kappa**(3/2)))) &
      * Ebin_keV * (1+Ebin_keV/(kappa*E0_keV))**(-1-1*kappa)
  case (2) ! multiple Maxwellians
    if (bimax_frac <= 0) then
      error stop 'ERROR:ionize_fang:fang2010_spectrum: bimax_frac must be greater than 0'
    end if
    phi_keV = (Q0_keV/(2*E0_keV**3)) * (Ebin_keV/bimax_frac**2) * exp(-1*Ebin_keV/(bimax_frac*E0_keV)) ! secondary Maxwellian
    phi_keV = phi_keV + (Q0_keV/(2*E0_keV**3)) * Ebin_keV * exp(-1*Ebin_keV/E0_keV) ! add primary Maxwellian
    phi_keV = phi_keV/2 ! normalize
  case (3) ! accelerated Maxwellian
    if (E0_char_keV < 0) then
      error stop 'ERROR:ionize_fang:fang2010_spectrum: W0_char must be greater than 0'
    end if
    if (Ebin_keV<E0_keV) then ! no electrons with energies less than U0
      phi_keV = 0
    else ! E0_kev = acc. potential, E0_char_keV = thermal energy
      phi_keV = (Q0_keV/(E0_char_keV**2+(E0_char_keV+E0_keV)**2)) * (Ebin_keV/E0_char_keV) * exp(-1*(Ebin_keV-E0_keV)/E0_char_keV)
    end if
  case (4) ! Data used from Evans, D. S. (1974)
    phi_keV = (Q0_keV/4.54615e9_wp)*phi_Evans_keV(k) ! normalizes to Q0 map
  case default
    error stop 'ERROR:ionize_fang:fang2010_spectrum: unknown diff_num_flux'
  end select

  !! differential precipitating energy flux
  dQ0_keV = Ebin_keV * phi_keV * dEbin_keV ! keV/s/cm^2

  !! Equation (3)
  !! total ionization rate "qtot" [cm^-3 s^-1]
  qtot = qtot + f * dQ0_keV / deps / H_cm
end do

end function fang2010_spectrum


elemental real(wp) function fang2010(Q0_keV, Emono_keV, Tn, massden_gcm3, meanmass_g, g_ms2) result(qtot)
!! isotropically precipitating monoenergetic (100 eV to 1 MeV) electrons
!! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GL045406

real(wp), intent(in) :: Q0_keV, Emono_keV, Tn, massden_gcm3, meanmass_g, g_ms2

real(wp) :: y, H_cm, f
integer :: i, j
real(wp), dimension(8) :: C
character(12) :: E0_str

real(wp), parameter :: P(8,4) = reshape( &
[1.24616_wp,     1.45903_wp,    -2.42269e-1_wp,  5.95459e-2_wp, &
 2.23976_wp,    -4.22918e-7_wp,  1.36458e-2_wp,  2.53332e-3_wp, &
 1.41754_wp,     1.44597e-1_wp,  1.70433e-2_wp,  6.39717e-4_wp, &
 2.48775e-1_wp, -1.50890e-1_wp,  6.30894e-9_wp,  1.23707e-3_wp, &
-4.65119e-1_wp, -1.05081e-1_wp, -8.95701e-2_wp,  1.22450e-2_wp, &
 3.86019e-1_wp,  1.75430e-3_wp, -7.42960e-4_wp,  4.60881e-4_wp, &
-6.45454e-1_wp,  8.49555e-4_wp, -4.28581e-2_wp, -2.99302e-3_wp, &
 9.48930e-1_wp,  1.97385e-1_wp, -2.50660e-3_wp, -2.06938e-3_wp], shape(P), order=[2,1])

if (Emono_keV < 0.1_wp .or. Emono_keV > 1000) then
  write(E0_str,'(F12.4)') Emono_keV
  error stop 'ionize_fang:fang2010: valid E0 range from 100 eV .. 1 MeV: E0 (keV) : ' // E0_str
endif

!! scale height
!! Equation (2)
H_cm = 100 * kb * Tn / (meanmass_g/1000) / abs(g_ms2)

!! normalized atmospheric column mass
!! Equation (1)
y = 2/Emono_keV * (massden_gcm3 * H_cm / 6e-6_wp)**0.7_wp

!! Equation (5)
C = 0
do i=1,size(P,1)
  do j=1,size(P,2)
    C(i) = C(i) + P(i,j) * log(Emono_keV)**(j-1)
  end do
end do
C = exp(C)

!! Equation (4)
!! Energy deposition "f"
f = C(1)*y**C(2)*exp(-1*C(3)*y**C(4)) + C(5)*y**C(6) * exp(-1*C(7)*y**C(8))

!! Equation (3)
!! total ionization rate "qtot" [cm^-3 s^-1]
qtot = f * Q0_keV / deps / H_cm

end function fang2010


elemental real(wp) function fang2008(Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g_ms2) result(qtot)

!! COMPUTE IONIZATION RATES PER THE FANG 2008 SEMI-EMPIRICAL METHOD.
!! https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384

!! Total Ionization Rate by Precipitating Electrons With a
!! Maxwellian Energy and Isotropic Pitch Angle Distribution

!! valid range of E0: 100 eV to 1 MeV

real(wp), intent(in) :: Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g_ms2
!! Q0: [keV cm^-2 s^-1]
!! massden_gcm3: [g cm^-3]
!! meanmass_g: [g]
!! g_ms2: [m s^-2]

real(wp) :: y, H_cm, f
integer :: i, j

character(12) :: E0_str

real(wp), dimension(8) :: C

real(wp), parameter :: P(8,4) = reshape( &
[3.49979e-1_wp, -6.18200e-2_wp, -4.08124e-2_wp,  1.65414e-2_wp, &
 5.85425e-1_wp, -5.00793e-2_wp,  5.69309e-2_wp, -4.02491e-3_wp, &
 1.69692e-1_wp, -2.58981e-2_wp,  1.96822e-2_wp,  1.20505e-3_wp, &
-1.22271e-1_wp, -1.15532e-2_wp,  5.37951e-6_wp,  1.20189e-3_wp, &
 1.57018_wp,     2.87896e-1_wp, -4.14857e-1_wp,  5.18158e-2_wp, &
 8.83195e-1_wp,  4.31402e-2_wp, -8.33599e-2_wp,  1.02515e-2_wp, &
 1.90953_wp,    -4.74704e-2_wp, -1.80200e-1_wp,  2.46652e-2_wp, &
-1.29566_wp,    -2.10952e-1_wp,  2.73106e-1_wp, -2.92752e-2_wp], shape(P), order=[2,1])

if (E0_keV < 0.1_wp .or. E0_keV > 1000) then
  write(E0_str,'(F12.4)') E0_keV
  error stop 'ionize_fang:fang2008: valid E0 range from 100 eV .. 1 MeV: E0 (keV) : ' // E0_str
endif

!! Equation (3)
!! scale height
!! kb [m^2 kg s^-2 K^-1]
H_cm = 100 * kb * Tn / (meanmass_g/1000) / abs(g_ms2)
!! H_cm [centimeters]

!! Equation (4)
y = 1 / E0_keV*(massden_gcm3 * H_cm / 4e-6_wp)**(0.606_wp)


!! Ci COEFFS and SHAPE FUNCTION
!! Equation (7)

C = 0
do i=1,size(P,1)
  do j=1,size(P,2)
    C(i) = C(i) + P(i,j) * log(E0_keV)**(j-1)
  end do
end do
C = exp(C)

!! Equation (6) energy deposition "f"
f = C(1)*y**C(2)*exp(-1*C(3)*y**C(4))+C(5)*y**C(6)*exp(-1*C(7)*y**C(8))

!! Equation (2) total electron impact ionization rate
qtot = Q0_keV / 2._wp / deps / H_cm * f
!! [cm^-3 s^-1]

end function fang2008


elemental real(wp) function gravity_accel(alt_km)
!! computes gravitational acceleration normal to Earth, where up is positive acceleration
real(wp), intent(in) :: alt_km
gravity_accel = -1 * 6.67408e-11_wp * 5.9722e24_wp / (6371.e3_wp + alt_km*1000)**2
!! [m s^-2]
end function gravity_accel

end module ionize_fang
