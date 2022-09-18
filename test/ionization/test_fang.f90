program test_fang
!! Reproduces data of:
!! * Figure 3 in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JA013384
!! * Figure 2 in Fang 2010

use, intrinsic:: iso_fortran_env, only: real32, real64
use assert, only : assert_isclose
use phys_consts, only: wp
use ionrate, only: ionization_fang2008, ionization_fang2010
use msis_interface, only : msisinit

implicit none (type, external)

integer :: i, argc
real(wp) :: Q0_erg, f107, f107a, Ap, glat, glon, UTsec, altrange(3), a
real(wp), allocatable :: alt_km(:), E0_keV(:), Qtot08(:,:), Qtot10(:,:)
integer :: doy, msis_version
real(wp) :: massden_gcm3, meanmass_g
real(wp), dimension(2) :: Q08ref, Q10ref

!! default parameters if no arguments given on CLI

Q0_erg = 1
E0_keV = [0.1, 1., 10., 100., 1000.]
altrange = [20, 400, 2]
f107=50
f107a=50
Ap=5
glat=60
glon=0
doy=1
UTsec=0
msis_version=0

Q08ref = 0
Q10ref = 0
!! dummy value

argc = command_argument_count()

if(argc > 0) call parse_arg(1, msis_version)

if (argc < 2) then
  select case (msis_version)
  case (0)
    Q08ref = [2214.052_wp, 9579.046_wp]
    Q10ref = [1192.002_wp, 778.655_wp]
  case (21)
    Q08ref = [2214.052_wp, 9579.046_wp]
    Q10ref = [1192.002_wp, 778.655_wp]
  end select
else
  call parse_arg(2, Q0_erg)

  deallocate(E0_keV)
  allocate(E0_keV(1))
  call parse_arg(3, E0_keV(1))

  do i = 4,6
    call parse_arg(i, altrange(i-3))
  enddo

  call parse_arg(7, f107)
  call parse_arg(8, f107a)
  call parse_arg(9, Ap)
  call parse_arg(10, glat)
  call parse_arg(11, glon)
  call parse_arg(12, doy)
  call parse_arg(13, UTsec)
endif

!! main program

if(msis_version > 0) then
  print '(A,f3.1)', "test_fang: initialized MSIS ", real(msis_version) / 10
  call msisinit()
endif

i = int((altrange(2)-altrange(1)) / altrange(3)) + 1
allocate(alt_km(i))
alt_km(1) = altrange(1)
do i = 2,size(alt_km)
  alt_km(i) = alt_km(i-1) + altrange(3)
enddo

allocate(Qtot08(size(alt_km), size(E0_keV)), Qtot10(size(alt_km), size(E0_keV)))


do i = 1, size(alt_km)
  Qtot08(i, :) = ionization_fang2008(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec, &
    msis_version=msis_version)
  Qtot10(i, :) = ionization_fang2010(Q0_erg, E0_keV, alt_km(i), f107, f107a, Ap, glat, glon, doy, UTsec, &
    msis_version=msis_version)
enddo

print '(A,25F10.1)', 'alt[km]/E0[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot08(i,:)
enddo

print '(/,A,25F10.1)', 'alt[km]/Emono[keV]', E0_keV
do i = 1,size(alt_km)
  print '(F7.1,5F15.3)', alt_km(i), Qtot10(i,:)
enddo

if (argc < 2) call compare_out(Qtot08, Qtot10, Q08ref, Q10ref)


contains


subroutine compare_out(Q08, Q10, Q08ref, Q10ref)
!! compare against arbitrary reference point data for known inputs
real(wp), intent(in), dimension(:,:) :: Q08, Q10
real(wp), intent(in), dimension(:) :: Q08ref, Q10ref

call assert_isclose(Q08(90, 1), Q08ref(1), atol=0.001_wp, err_msg="E0: 100eV")
call assert_isclose(Q08(18, 5), Q08ref(2), atol=0.001_wp, err_msg="E0: 1MeV")

call assert_isclose(Q10(90, 1), Q10ref(1), atol=0.001_wp, err_msg="Emono: 100eV")
call assert_isclose(Q10(18, 5), Q10ref(2), atol=0.001_wp, rtol=0.001_wp, err_msg="Emono: 1MeV")

end subroutine compare_out


subroutine parse_arg(iarg, val)

integer, intent(in) :: iarg
class(*), intent(out) :: val

integer :: ierr
character(80) :: argv

call get_command_argument(iarg, argv, status=ierr)
if(ierr /= 0) error stop "ERROR: parse_args: could not get argument"

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
  error stop 'ERROR: parse_args: unknown type to parse'
end select

end subroutine parse_arg


end program
