program test_fang
!! Reproduces data of:
!! * Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
!! * Figure 2 in Fang 2010

use, intrinsic:: iso_fortran_env, only: real32, real64
use assert, only : assert_isclose
use phys_consts, only: wp
use ionrate, only: ionization_fang2008, ionization_fang2010

implicit none (type, external)

integer :: i, argc
real(wp) :: Q0_erg, f107, f107a, Ap, glat, glon, UTsec, altrange(3)
real(wp), allocatable :: alt_km(:), E0_keV(:), Qtot08(:,:), Qtot10(:,:)
integer :: doy

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

i = int((altrange(2)-altrange(1)) / altrange(3)) + 1
allocate(alt_km(i))
alt_km(1) = altrange(1)
do i = 2,size(alt_km)
  alt_km(i) = alt_km(i-1) + altrange(3)
enddo

allocate(Qtot08(size(alt_km), size(E0_keV)), Qtot10(size(alt_km), size(E0_keV)))


do i = 1, size(alt_km)
  Qtot08(i, :) = ionization_fang2008(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec, msis_version=0)
  Qtot10(i, :) = ionization_fang2010(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec, msis_version=0)
enddo

print '(A,25F10.1)', 'alt[km]/E0[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot08(i,:)
enddo

print '(/,A,25F10.1)', 'alt[km]/Emono[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot10(i,:)
enddo

call assert_isclose(Qtot08(90, 1), 2214.052_wp, atol=0.001_wp, err_msg="E0: 100eV")
call assert_isclose(Qtot08(18, 5), 9579.046_wp, atol=0.001_wp, err_msg="E0: 1MeV")

call assert_isclose(Qtot10(90, 1), 1192.002_wp, atol=0.001_wp, err_msg="Emono: 100eV")
call assert_isclose(Qtot10(18, 5), 778.655_wp, atol=0.001_wp, rtol=0.001_wp, err_msg="Emono: 1MeV")

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
