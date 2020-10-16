Program test_fang
!! Reproduces data of:
!! * Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
!! * Figure 2 in Fang 2010

use, intrinsic:: iso_fortran_env, only: real32, real64
use phys_consts, only: wp
use ionrate, only: ionization_fang2008, ionization_fang2010
use ionize_fang, only : maxwellian_dnf, erg2kev

implicit none (type, external)

integer :: i, argc
real(wp) :: Q0_erg, f107, f107a, Ap, glat, glon, UTsec, altrange(3), a
real(wp), allocatable :: alt_km(:), E0_keV(:), Qtot08(:,:), Qtot10(:,:), Qtot10_fig4(:)
integer :: doy
real(wp) :: massden_gcm3, meanmass_g

!> Energy bins [keV]
!! numpy.logspace(-1, 3, 200)
real(wp), parameter :: Ebin(*) = &
[1.00000000e-01, 1.04737090e-01, 1.09698580e-01, 1.14895100e-01, &
1.20337784e-01, 1.26038293e-01, 1.32008840e-01, 1.38262217e-01, &
1.44811823e-01, 1.51671689e-01, 1.58856513e-01, 1.66381689e-01, &
1.74263339e-01, 1.82518349e-01, 1.91164408e-01, 2.00220037e-01, &
2.09704640e-01, 2.19638537e-01, 2.30043012e-01, 2.40940356e-01, &
2.52353917e-01, 2.64308149e-01, 2.76828663e-01, 2.89942285e-01, &
3.03677112e-01, 3.18062569e-01, 3.33129479e-01, 3.48910121e-01, &
3.65438307e-01, 3.82749448e-01, 4.00880633e-01, 4.19870708e-01, &
4.39760361e-01, 4.60592204e-01, 4.82410870e-01, 5.05263107e-01, &
5.29197874e-01, 5.54266452e-01, 5.80522552e-01, 6.08022426e-01, &
6.36824994e-01, 6.66991966e-01, 6.98587975e-01, 7.31680714e-01, &
7.66341087e-01, 8.02643352e-01, 8.40665289e-01, 8.80488358e-01, &
9.22197882e-01, 9.65883224e-01, 1.01163798e+00, 1.05956018e+00, &
1.10975250e+00, 1.16232247e+00, 1.21738273e+00, 1.27505124e+00, &
1.33545156e+00, 1.39871310e+00, 1.46497140e+00, 1.53436841e+00, &
1.60705282e+00, 1.68318035e+00, 1.76291412e+00, 1.84642494e+00, &
1.93389175e+00, 2.02550194e+00, 2.12145178e+00, 2.22194686e+00, &
2.32720248e+00, 2.43744415e+00, 2.55290807e+00, 2.67384162e+00, &
2.80050389e+00, 2.93316628e+00, 3.07211300e+00, 3.21764175e+00, &
3.37006433e+00, 3.52970730e+00, 3.69691271e+00, 3.87203878e+00, &
4.05546074e+00, 4.24757155e+00, 4.44878283e+00, 4.65952567e+00, &
5.87278661e+00, 6.15098579e+00, 6.44236351e+00, 6.74754405e+00, &
7.06718127e+00, 7.40196000e+00, 7.75259749e+00, 8.11984499e+00, &
8.50448934e+00, 8.90735464e+00, 9.32930403e+00, 9.77124154e+00, &
1.02341140e+01, 1.07189132e+01, 1.12266777e+01, 1.17584955e+01, &
1.23155060e+01, 1.28989026e+01, 1.35099352e+01, 1.41499130e+01, &
1.48202071e+01, 1.55222536e+01, 1.62575567e+01, 1.70276917e+01, &
1.78343088e+01, 1.86791360e+01, 1.95639834e+01, 2.04907469e+01, &
2.14614120e+01, 2.24780583e+01, 2.35428641e+01, 2.46581108e+01, &
2.58261876e+01, 2.70495973e+01, 2.83309610e+01, 2.96730241e+01, &
3.10786619e+01, 3.25508860e+01, 3.40928507e+01, 3.57078596e+01, &
3.73993730e+01, 3.91710149e+01, 4.10265811e+01, 4.29700470e+01, &
4.50055768e+01, 4.71375313e+01, 4.93704785e+01, 5.17092024e+01, &
5.41587138e+01, 5.67242607e+01, 5.94113398e+01, 6.22257084e+01, &
6.51733960e+01, 6.82607183e+01, 7.14942899e+01, 7.48810386e+01, &
7.84282206e+01, 8.21434358e+01, 8.60346442e+01, 9.01101825e+01, &
9.43787828e+01, 9.88495905e+01, 1.03532184e+02, 1.08436597e+02, &
1.13573336e+02, 1.18953407e+02, 1.24588336e+02, 1.30490198e+02, &
1.36671636e+02, 1.43145894e+02, 1.49926843e+02, 1.57029012e+02, &
1.64467618e+02, 1.72258597e+02, 1.80418641e+02, 1.88965234e+02, &
1.97916687e+02, 2.07292178e+02, 2.17111795e+02, 2.27396575e+02, &
2.38168555e+02, 2.49450814e+02, 2.61267523e+02, 2.73644000e+02, &
2.86606762e+02, 3.00183581e+02, 3.14403547e+02, 3.29297126e+02, &
3.44896226e+02, 3.61234270e+02, 3.78346262e+02, 3.96268864e+02, &
4.15040476e+02, 4.34701316e+02, 4.55293507e+02, 4.76861170e+02, &
4.99450512e+02, 5.23109931e+02, 5.47890118e+02, 5.73844165e+02, &
6.01027678e+02, 6.29498899e+02, 6.59318827e+02, 6.90551352e+02, &
7.23263390e+02, 7.57525026e+02, 7.93409667e+02, 8.30994195e+02, &
8.70359136e+02, 9.11588830e+02, 9.54771611e+02, 1.00000000e+03]

real(wp) :: double_maxwellian_dnf(size(Ebin))

argc = command_argument_count()
if (argc < 1) then
  Q0_erg = 1
else
  call parse_arg(1, Q0_erg)
endif
if (argc < 2) then
  E0_keV = [0.1, 1., 10., 100., 1000.]
else
  allocate(E0_keV(1))
  call parse_arg(2, E0_keV(1))
endif
if (argc<5) then
  altrange = [20, 400, 2]
else
  do i = 3,5
    call parse_arg(i, altrange(i-2))
  enddo
endif
if (argc<12) then
  f107=50
  f107a=50
  Ap=5
  glat=60
  glon=0
  doy=1
  UTsec=0
else
  call parse_arg(6, f107)
  call parse_arg(7, f107a)
  call parse_arg(8, Ap)
  call parse_arg(9, glat)
  call parse_arg(10, glon)
  call parse_arg(11, doy)
  call parse_arg(12, UTsec)
endif

!> altitude grid

i = int((altrange(2)-altrange(1)) / altrange(3)) + 1
allocate(alt_km(i))
alt_km(1) = altrange(1)
do i = 2,size(alt_km)
  alt_km(i) = alt_km(i-1) + altrange(3)
enddo

!> Figure 2 and 3

allocate(Qtot08(size(alt_km), size(E0_keV)), Qtot10(size(alt_km), size(E0_keV)))

do i = 1, size(alt_km)
  Qtot08(i, :) = ionization_fang2008(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec)
  Qtot10(i, :) = ionization_fang2010(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec)
enddo

print '(A,25F10.1)', 'alt[km]/E0[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot08(i,:)
enddo

print '(/,A,25F10.1)', 'alt[km]/Emono[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot10(i,:)
enddo

!> Figure 4

allocate(Qtot10_fig4(size(alt_km)))

double_maxwellian_dnf = (maxwellian_dnf(Q0_keV=1._wp * erg2kev, E0_keV=5._wp, E=Ebin) + &
                         maxwellian_dnf(Q0_keV=0.5_wp*erg2kev, E0_keV=50._wp, E=Ebin)) / erg2kev

do i = 1, size(alt_km)
  Qtot10_fig4(i) = sum(ionization_fang2010(double_maxwellian_dnf, Ebin, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec))
enddo

print '(/,A,500F10.3)', 'Ebin[keV]', Ebin
print '(A,500F10.3)', 'Qflux [keV^-1 cm^-2 s^-1]', double_maxwellian_dnf
print *, 'alt[km], IonizationRate'
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot10_fig4(i)
enddo

contains

subroutine parse_arg(iarg, val)
  integer, intent(in) :: iarg
  class(*), intent(out) :: val

  character(80) :: argv

  call get_command_argument(iarg, argv)

  select type(val)
    type is (integer)
      read(argv,*) val
    type is (real(real64))
      read(argv,*) val
    type is (real(real32))
      read(argv,*) val
    type is (character(*))
      read(argv,*) val
    class default
      error stop 'unknown type to parse'
  end select

end subroutine parse_arg


end program
