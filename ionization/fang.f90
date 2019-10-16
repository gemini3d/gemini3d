module ionize_fang

use phys_consts, only: wp, kb
implicit none
private

real(wp), parameter :: deps = 0.035_wp
!! keV, kinetic energy lost per ion-electron pair produced

public :: fang2008, fang2010

contains

elemental real(wp) function fang2010(Q0, Emono_keV, Tn, massden_gcm3, meanmass, g1) result(Qtot)

real(wp), intent(in) :: Emono_keV, Q0, Tn, massden_gcm3, meanmass, g1

real(wp) :: y, H_cm, f
integer :: i, j
real(wp), dimension(8) :: C


real(wp), parameter :: P(8,4) = reshape( &
[1.24616_wp,     1.45903_wp,    -2.42269e-1_wp,  5.95459e-2_wp, &
 2.23976_wp,    -4.22918e-7_wp,  1.36458e-2_wp,  2.53332e-3_wp, &
 1.41754_wp,     1.44597e-1_wp,  1.70433e-2_wp,  6.39717e-4_wp, &
 2.48775e-1_wp, -1.50890e-1_wp,  6.30894e-9_wp,  1.23707e-3_wp, &
-4.65119e-1_wp, -1.05081e-1_wp, -8.95701e-2_wp,  1.22450e-2_wp, &
 3.86019e-1_wp,  1.75430e-3_wp, -7.42960e-4_wp,  4.60881e-4_wp, &
-6.45454e-1_wp,  8.49555e-4_wp, -4.28581e-2_wp, -2.99302e-3_wp, &
 9.48930e-1_wp,  1.97385e-1_wp, -2.50660e-3_wp, -2.06938e-3_wp], shape(P), order=[2,1])


!! scale height
!! Equation (2)
H_cm = kb * Tn / meanmass / abs(g1)

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
Qtot = f * Q0 / deps / H_cm

end function fang2010


impure elemental real(wp) function fang2008(Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g1_ms2) result(Qtot)

!! COMPUTE IONIZATION RATES PER THE FANG 2008 SEMI-EMPIRICAL METHOD.
!! https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384

!! Total Ionization Rate by Precipitating Electrons With a
!! Maxwellian Energy and Isotropic Pitch Angle Distribution

real(wp), intent(in) :: Q0_keV, E0_keV, Tn, massden_gcm3, meanmass_g, g1_ms2
!! Q0: [keV cm^-2 s^-1]
!! massden_gcm3: [g cm^-3]
!! meanmass_g: [g]
!! g1_ms2: [m s^-2]

real(wp) :: y, H_cm, f
integer :: i, j

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

!! Equation (3)
!! scale height
!! kb [m^2 kg s^-2 K^-1]
H_cm = 100 * kb * Tn / (meanmass_g/1000) / abs(g1_ms2)
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
Qtot = Q0_keV / 2._wp / deps / H_cm * f
!! [cm^-3 s^-1]

!print *, 'massden, meanmass:',massden_gcm3, meanmass_g
!print *,'y',y

end function fang2008

end module ionize_fang