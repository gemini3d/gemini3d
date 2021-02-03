module assert

use, intrinsic:: iso_fortran_env, only: stderr=>error_unit, real64, real32
use, intrinsic:: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan

implicit none (type, external)

private

public :: isclose, assert_allclose

contains

elemental logical function isclose(actual, desired, rtol, atol, equal_nan)
!! ## inputs
!!
!! * actual: value "measured"
!! * desired: value "wanted"
!! * rtol: relative tolerance
!! * atol: absolute tolerance
!! * equal_nan: consider NaN to be equal?
!!
!!  rtol overrides atol when both are specified
!!
!! https://www.python.org/dev/peps/pep-0485/#proposed-implementation
!! https://github.com/PythonCHB/close_pep/blob/master/is_close.py

class(*), intent(in) :: actual, desired
class(*), intent(in), optional :: rtol, atol
logical, intent(in), optional :: equal_nan

real(real64) :: r,a, act, des
logical :: n

isclose = .false. !< ensure it's defined

!> INSTEAD OF merge(), since non present values aren't defined.
r = 1e-6
a = 1e-12
n = .false.

select type (actual)
type is (real(real32))
  act = real(actual, real32)
type is (real(real64))
  act = actual
class default
  error stop "assert: actual must be real32 or real64"
end select

select type (desired)
type is (real(real32))
  des = real(desired, real32)
type is (real(real64))
  des = desired
class default
  error stop "assert: desired must be real32 or real64"
end select

if (present(rtol)) then
  select type (rtol)
  type is (real(real64))
    r = rtol
  type is (real(real32))
    r = real(rtol, real64)
  class default
    error stop "assert: rtol needs real32 or real64"
  end select
endif

if (present(atol)) then
  select type (atol)
  type is (real(real64))
    a = atol
  type is (real(real32))
    a = real(atol, real64)
  class default
    error stop "assert: atol needs real32 or real64"
  end select
endif

if (present(equal_nan)) n = equal_nan

!print*,r,a,n,act,des

!> sanity check
if ((r < 0).or.(a < 0)) error stop 'improper tolerances specified'
!> simplest case -- too unlikely, especially for arrays?
!isclose = (act == des)
!if (isclose) return

!> equal nan
if (n) then ! fortran is NOT short circuit logic in general
  isclose = (ieee_is_nan(act) .and. ieee_is_nan(des))
  if (isclose) return
endif

!> Inf /= -Inf, unequal NaN
if (.not.ieee_is_finite(act) .or. .not.ieee_is_finite(des)) return

!> floating point closeness check
isclose = abs(act-des) <= max(r * max(abs(act), abs(des)), a)

end function isclose


impure elemental subroutine assert_allclose(actual, desired, rtol, atol, equal_nan, err_msg)

!! ## inputs
!!
!! *  actual: value "measured"
!! *  desired: value "wanted"
!! *  rtol: relative tolerance
!! *  atol: absolute tolerance
!! *  equal_nan: consider NaN to be equal?
!! *  err_msg: message to print on mismatch
!!
!! rtol overrides atol when both are specified

class(*), intent(in) :: actual, desired
class(*), intent(in), optional :: rtol, atol
logical, intent(in), optional :: equal_nan
character(*), intent(in), optional :: err_msg
character(:), allocatable :: emsg

real(real64) :: act, des

select type (actual)
type is (real(real32))
  act = real(actual, real32)
type is (real(real64))
  act = actual
class default
  error stop "assert: actual must be real32 or real64"
end select

select type (desired)
type is (real(real32))
  des = real(desired, real32)
type is (real(real64))
  des = desired
class default
  error stop "assert: desired must be real32 or real64"
end select

if (present(err_msg)) then
  emsg = err_msg
else
  emsg = ''
endif

if (.not.isclose(act, des, rtol,atol,equal_nan)) then
  write(stderr,*) emsg // ': actual',act,'desired',des
  error stop
endif

end subroutine assert_allclose

end module assert
