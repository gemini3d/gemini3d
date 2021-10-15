module msis_interface
!! this module allows selecting MSISE00 or MSIS 2.0
!! it is a thin abstraction of the MSIS routines
!! MSISE00 is Fortran 66 style, while MSIS 2.0 is Fortran 90 style
!!
!! We assume MSISE00 is always available, which MSIS 2.0 might not be available.

use msis_calc, only : msiscalc
use msis_init, only : msisinit
use, intrinsic :: iso_fortran_env, only : real32, real64, int32
implicit none (type, external)

interface msis_gtd7
  module procedure msis_gtd7_r32, msis_gtd7_r64
end interface msis_gtd7

interface msis_gtd8
  module procedure msis_gtd8_r64, msis_gtd8_r32
end interface msis_gtd8

private
public :: msis_gtd7, msis_gtd8, msisinit

contains

subroutine msis_gtd7_r32(doy, UTsec, alt_km,  glat, glon, f107a, f107, Ap7, d, T, use_meters, sw25)

use msise00_gemini, only : meters, gtd7, tselec

class(*), intent(in) :: doy
real(real32), intent(in) :: UTsec, alt_km, glat, glon, f107, f107a, Ap7(7)
real(real32), intent(out) :: d(9),T(2)
logical, intent(in) :: use_meters
real(real32), intent(in), optional :: sw25(25)

real(real32) :: stl, sw(25)
integer :: dayOfYear

select type (doy)
  type is (integer(int32))
    dayOfYear = doy
  type is (real(real64))
    dayOfYear = int(doy)
  type is (real(real32))
    dayOfYear = int(doy)
  class default
    error stop 'msis_gtd8: doy must be real or integer'
end select

!> input validation
if (dayOfYear < 1 .or. dayOfYear > 366) error stop "valid dayOfYear range is 1..366"
if (UTsec < 0 .or. UTsec > 86400) error stop "valid UTsec range is 0..86400"

!> input conversion
stl = UTsec/3600 + glon/15

! print '(A,I0.3,14F8.1)', 'TRACE:MSIS00_bit32: inputs: ',dayOfYear, UTsec, alt_km, glat, glon, stl, f107a, f107, Ap7

call meters(use_meters)

sw = 1
if (present(sw25)) sw = sw25
! print '(A,25F4.1)','TRACE:MSIS00_bit32: sw25: ',sw
call tselec(sw)

call gtd7(dayOfYear, UTsec, alt_km, glat, glon, stl, f107a, f107, Ap7, 48, d, T)

! print *,'TRACE:MSIS00_bit32: d,T:', d,T

end subroutine msis_gtd7_r32


subroutine msis_gtd7_r64(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, d, T, use_meters, sw25)
!! adds casting to/from real32
use msise00_gemini, only : meters, gtd7, tselec

class(*), intent(in) :: doy
real(real64), intent(in) :: UTsec, alt_km, glat, glon, f107, f107a, Ap7(7)
real(real64), intent(out) :: d(9),T(2)
logical, intent(in) :: use_meters
real(real64), intent(in), optional :: sw25(25)

real(real32) :: sw(25), stl, d32(9), T32(2)
integer :: dayOfYear

select type (doy)
  type is (integer(int32))
    dayOfYear = doy
  type is (real(real64))
    dayOfYear = int(doy)
  type is (real(real32))
    dayOfYear = int(doy)
  class default
    error stop 'msis_gtd8: doy must be real or integer'
end select

!> input validation
if (dayOfYear < 1 .or. dayOfYear > 366) error stop "valid dayOfYear range is 1..366"
if (UTsec < 0 .or. UTsec > 86400) error stop "valid UTsec range is 0..86400"

!> input conversion
stl = real(UTsec/3600 + glon/15, real32)

call meters(use_meters)

sw = 1
if (present(sw25)) sw = real(sw25, real32)
call tselec(sw)

call gtd7(dayOfYear, real(UTsec, real32), real(alt_km, real32), &
  real(glat, real32), real(glon, real32), real(stl, real32), &
  real(f107a, real32), real(f107, real32), real(Ap7, real32), 48, &
  d32, T32)

d = real(d32, real64)
T = real(T32, real64)

end subroutine msis_gtd7_r64


subroutine msis_gtd8_r64(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, Dn, Tn)
!! translate MSIS 2.0 to MSISE00 gtd7-like
!! assume MSIS 2.0 is real32
class(*), intent(in) :: doy
real(real64), intent(in) :: UTsec, alt_km, glat, glon, f107a, f107, Ap7(7)
real(real64), intent(out) :: Dn(9), Tn(2)

real(real32) :: D(10), T(2), dayOfYear

select type (doy)
  type is (integer(int32))
    dayOfYear = real(doy, real32)
  type is (real(real64))
    dayOfYear = real(doy, real32)
  type is (real(real32))
    dayOfYear = doy
  class default
    error stop 'msis_gtd8: doy must be real or integer'
end select

! print *, "TRACE: MSIS 2.0 64 bit: inputs: doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7", &
!   dayOfYear, UTsec, alt_km, glat, glon, f107a, f107, Ap7(1)

!> input validation
if (dayOfYear < 1 .or. dayOfYear > 366) error stop "valid dayOfYear range is 1..366"
if (UTsec < 0 .or. UTsec > 86400) error stop "valid UTsec range is 0..86400"

call msiscalc(day=dayOfYear, UTsec=real(UTsec, real32), &
  z=real(alt_km, real32), lat=real(glat, real32), lon=real(glon, real32), &
  SfluxAvg=real(f107a, real32), Sflux=real(f107, real32), ap=real(Ap7, real32), &
  Tn=T(2), Tex=T(1), Dn=D)

!> translate to old gtd7 convention
Dn(1) = D(5)  !< He
Dn(2) = D(4)  !< O
Dn(3) = D(2)  !< N2
Dn(4) = D(3)  !< O2
Dn(5) = D(7)  !< Ar
Dn(6) = D(1)  !< Total mass density
Dn(7) = D(6)  !< H
Dn(8) = D(8)  !< N
Dn(9) = D(9)  !< Anomalous O
!! D(10) will be NO in future MSIS 2.x

Tn = real(T, real64)

end subroutine msis_gtd8_r64


subroutine msis_gtd8_r32(doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7, Dn, Tn)
!! translate MSIS 2.0 to MSISE00 gtd7-like
!! assume MSIS 2.0 is also real32
class(*), intent(in) :: doy
real(real32), intent(in) :: UTsec, alt_km, glat, glon, f107a, f107, Ap7(7)
real(real32), intent(out) :: Dn(9), Tn(2)

real(real32) :: D(10), dayOfYear

select type (doy)
  type is (integer(int32))
    dayOfYear = real(doy, real32)
  type is (real(real64))
    dayOfYear = real(doy, real32)
  type is (real(real32))
    dayOfYear = doy
  class default
    error stop 'msis_gtd8: doy must be real or integer'
end select

! print *, "TRACE: MSIS 2.0 64 bit: inputs: doy, UTsec, alt_km, glat, glon, f107a, f107, Ap7", &
!   dayOfYear, UTsec, alt_km, glat, glon, f107a, f107, Ap7(1)

!> input validation
if (dayOfYear < 1 .or. dayOfYear > 366) error stop "valid dayOfYear range is 1..366"
if (UTsec < 0 .or. UTsec > 86400) error stop "valid UTsec range is 0..86400"

call msiscalc(day=dayOfYear, UTsec=UTsec, &
  z=alt_km, lat=glat, lon=glon, &
  SfluxAvg=f107a, Sflux=f107, ap=Ap7, &
  Tn=Tn(2), Tex=Tn(1), Dn=D)

!> translate to old gtd7 convention
Dn(1) = D(5)  !< He
Dn(2) = D(4)  !< O
Dn(3) = D(2)  !< N2
Dn(4) = D(3)  !< O2
Dn(5) = D(7)  !< Ar
Dn(6) = D(1)  !< Total mass density
Dn(7) = D(6)  !< H
Dn(8) = D(8)  !< N
Dn(9) = D(9)  !< Anomalous O
!! D(10) will be NO in future MSIS 2.x

end subroutine msis_gtd8_r32

end module
