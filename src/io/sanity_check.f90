module sanity_check
!! check that variables are within realm of possible values,
!! at least that they're finite (not NaN or infinite)

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use phys_consts, only : wp

implicit none (type, external)
private
public :: check_finite_output, check_finite_plasma, check_finite_current, check_finite_mag

interface check_finite
module procedure check_1d, check_2d, check_3d, check_4d
end interface check_finite

contains


subroutine check_finite_output(time,id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

real(wp), intent(in) :: time
integer, intent(in) :: id
real(wp), intent(in), dimension(:,:,:,:) :: vs2, vs3, ns, vs1, Ts
real(wp), intent(in), dimension(:,:,:) :: Phi, J1, J2, J3

call check_finite(ns, 'ns', time, id)
call check_finite(vs1, 'vs1', time, id)
call check_finite(Ts, 'Ts', time, id)
call check_finite(J1, 'J1', time, id)
call check_finite(J2, 'J2', time, id)
call check_finite(J3, 'J3', time, id)
call check_finite(vs2, 'vs2', time, id)
call check_finite(vs3, 'vs3', time, id)
! call check_finite(Phi, 'Phi', time, id)

end subroutine check_finite_output


subroutine check_finite_plasma(ns, vs1, Ts)

real(wp), dimension(:,:,:,:), intent(in) :: ns
real(wp), dimension(:,:,:,:), intent(in) :: vs1, Ts

call check_finite(ns, 'ns')
call check_finite(vs1, 'vs1')
call check_finite(Ts, 'Ts')

if (any(ns < 0)) then
  write(stderr,'(A,ES15.3)') 'ERROR: negative number density minimum: ',minval(ns)
  error stop 'negative density'
endif
if (maxval(ns) < 1e6) then
  write(stderr,'(A,ES15.3)') 'ERROR: low number density maximum: ',maxval(ns)
  error stop 'unrealistically low maximum density'
endif
if (maxval(ns) > 1e16) then
  write(stderr,'(A,ES15.3)') 'ERROR: excessive number density maximum: ',maxval(ns)
  error stop 'unrealistically high maximum density'
endif

if (any(abs(vs1) > 1e7_wp)) then
  write(stderr,'(A,ES15.3)') 'ERROR: too fast velocity maximum: ',maxval(abs(vs1))
  error stop 'drift should not be realativistic'
endif

if (any(Ts < 0)) error stop 'negative temperature in Ts'
if (any(Ts > 100000)) error stop 'too hot Ts'
if (maxval(Ts) < 500) error stop 'too cold maximum Ts'

end subroutine check_finite_plasma


subroutine check_finite_mag(Br, Btheta, Bphi)

real(wp), dimension(:), intent(in) :: Br, Btheta, Bphi

call check_finite(Br, 'Br')
call check_finite(Btheta, 'Btheta')
call check_finite(Bphi, 'Bphi')

end subroutine check_finite_mag


subroutine check_finite_current(J1,J2,J3)

real(wp), intent(in), dimension(:,:,:) :: J1, J2, J3

call check_finite(J1, 'J1')
call check_finite(J2, 'J2')
call check_finite(J3, 'J3')

end subroutine check_finite_current


subroutine check_1d(var, name, time, id)

real(wp), intent(in) :: var(:)
character(*), intent(in) :: name
real(wp), intent(in), optional :: time
integer, intent(in), optional :: id

if (.not.all(ieee_is_finite(var))) call fail(name, time, id)

end subroutine check_1d


subroutine check_2d(var, name, time, id)

real(wp), intent(in) :: var(:,:)
character(*), intent(in) :: name
real(wp), intent(in), optional :: time
integer, intent(in), optional :: id

if (.not.all(ieee_is_finite(var))) call fail(name, time, id)

end subroutine check_2d


subroutine check_3d(var, name, time, id)

real(wp), intent(in) :: var(:,:,:)
character(*), intent(in) :: name
real(wp), intent(in), optional :: time
integer, intent(in), optional :: id

if (.not.all(ieee_is_finite(var))) call fail(name, time, id)

end subroutine check_3d


subroutine check_4d(var, name, time, id)

real(wp), intent(in) :: var(:,:,:,:)
character(*), intent(in) :: name
real(wp), intent(in), optional :: time
integer, intent(in), optional :: id

if (.not.all(ieee_is_finite(var))) call fail(name, time, id)

end subroutine check_4d


subroutine fail(name,time, id)

character(*), intent(in) :: name
real(wp), intent(in), optional :: time
integer, intent(in), optional :: id

write(stderr, '(A)', advance='no') 'sanity_check: non-finite value(s) in variable '// name

if (present(time) .and. present(id)) write(stderr,*) ' at time: ', time, ' on MPI worker # ', id

error stop

end subroutine fail


end module sanity_check