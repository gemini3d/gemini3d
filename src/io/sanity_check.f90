module sanity_check
!! check that variables are within realm of possible values,
!! at least that they're finite (not NaN or infinite)

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use phys_consts, only : wp
use errors, only : error_stop

implicit none (type, external)
private
public :: check_finite_output, check_finite_plasma, check_finite_current, check_finite_mag, check_finite_pertub

contains


pure subroutine ghost_bound(A, j1, k1, j2, k2, j3, k3, j4, k4)
!! "ig" is a priori the 2 ghost cells on each grid cell boundary, for MPI haloing.
real(wp), intent(in) :: A(..)
integer, intent(out) :: j1, k1, j2, k2, j3, k3
integer, intent(out), optional :: j4, k4

integer,parameter :: ig=2
integer :: r

r = rank(A)

if (r > 4 .or. r < 3) error stop "sanity_check:ghost_bound: only for rank 3,4 for now"
if (r == 4 .and. .not. (present(j4) .and. present(k4))) error stop "sanity_check:ghost_bound: 4d needs j4, k4"

j1 = lbound(A, 1) + ig
k1 = ubound(A, 1) - ig
j2 = lbound(A, 2) + ig
k2 = ubound(A, 2) - ig
j3 = lbound(A, 3) + ig
k3 = ubound(A, 3) - ig
if (r >= 4) then
  j4 = lbound(A, 4) + ig
  k4 = ubound(A, 4) - ig
endif

end subroutine ghost_bound


subroutine check_finite_output(out_dir, t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)
!! check outputs before proceeding to next time step
!! assumed to have GHOST CELLS, 2 at each end of each dimension
!! We purposely omit checking ghost cells, as they are for inter-worker communications
!! and so at the edges of the grid, the ghost cells contain random data.
!! Since the IEEE754 NaN allows for a wide range of values, it's probable that checking
!! large numbers of unassigned values will result in a non-finite value.

character(*), intent(in) :: out_dir
real(wp), intent(in) :: t_elapsed
integer, intent(in) :: worker_id
real(wp), intent(in), dimension(:,:,:,:) :: vs2, vs3, ns, vs1, Ts
real(wp), intent(in), dimension(:,:,:) :: Phi, J1, J2, J3

integer :: i1, k1, i2, k2, i3, k3, i4, k4
character(:), allocatable :: dump_filename
character(8) :: wid

write(wid, '(I0)') worker_id

dump_filename = out_dir // "/dump_nonfinite_output_worker_" // trim(wid) // ".h5"

call ghost_bound(ns, i1,k1, i2,k2, i3,k3, i4,k4)

if (.not.all(ieee_is_finite(ns(i1:k1, i2:k2, i3:k3, i4:k4)))) &
   call error_stop(dump_filename, 'output: non-finite Ns', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

call ghost_bound(vs1, i1,k1, i2,k2, i3,k3, i4,k4)

if (.not.all(ieee_is_finite(vs1(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'output: non-finite vs1', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

call ghost_bound(vs2, i1,k1, i2,k2, i3,k3, i4,k4)

if (.not.all(ieee_is_finite(vs2(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'output: non-finite vs2', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

call ghost_bound(vs3, i1,k1, i2,k2, i3,k3, i4,k4)

if (.not.all(ieee_is_finite(vs3(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'output: non-finite vs3', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

if (.not.all(ieee_is_finite(Ts(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'output: non-finite Ts', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

call ghost_bound(J1, i1,k1, i2,k2, i3,k3)

if (.not.all(ieee_is_finite(J1(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'output: non-finite J1', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

if (.not.all(ieee_is_finite(J2(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'output: non-finite J2', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

if (.not.all(ieee_is_finite(J3(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'output: non-finite J3', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

if (.not.all(ieee_is_finite(Phi(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'output: non-finite Phi', t_elapsed, worker_id, vs2,vs3,ns,vs1,Ts,Phi,J1,J2,J3)

end subroutine check_finite_output


subroutine check_finite_plasma(out_dir, ns, vs1, Ts)
!! check input data--garbage in, garbage out...

character(*), intent(in) :: out_dir
real(wp), dimension(:,:,:,:), intent(in) :: ns
real(wp), dimension(:,:,:,:), intent(in) :: vs1, Ts

integer :: i1, k1, i2, k2, i3, k3, i4, k4

character(:), allocatable :: dump_filename

dump_filename = out_dir // "/dump_nonfinite_plasma.h5"

call ghost_bound(ns, i1,k1, i2,k2, i3,k3, i4,k4)

if (.not.all(ieee_is_finite(ns(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'input:plasma: non-finite Ns', ns, vs1, Ts)

if (.not.all(ieee_is_finite(vs1(i1:k1, i2:k2, i3:k3, i4:k4)))) &
   call error_stop(dump_filename, 'input:plasma: non-finite vs1', ns, vs1, Ts)

if (.not.all(ieee_is_finite(Ts(i1:k1, i2:k2, i3:k3, i4:k4)))) &
  call error_stop(dump_filename, 'input:plasma: non-finite Ts', ns, vs1, Ts)

if (any(ns(i1:k1, i2:k2, i3:k3, i4:k4) < 0)) &
  call error_stop(dump_filename, 'input:plasma: negative density Ns', ns, vs1, Ts)

if (maxval(ns(i1:k1, i2:k2, i3:k3, i4:k4)) < 1e3) &
  call error_stop(dump_filename, 'input:plasma: too low maximum density', ns, vs1, Ts)

if (maxval(ns(i1:k1, i2:k2, i3:k3, i4:k4)) > 1e16) &
  call error_stop(dump_filename, 'input:plasma: too high maximum density', ns, vs1, Ts)

if (any(abs(vs1(i1:k1, i2:k2, i3:k3, i4:k4)) > 1e7_wp)) &
  call error_stop (dump_filename, 'input:plasma: drift realativistic', ns, vs1, Ts)

if (any(Ts(i1:k1, i2:k2, i3:k3, i4:k4) < 0)) &
  call error_stop (dump_filename, 'input:plasma: negative temperature in Ts', ns, vs1, Ts)

if (any(Ts(i1:k1, i2:k2, i3:k3, i4:k4) > 500000)) &
  call error_stop (dump_filename, 'input:plasma: too hot Ts', ns, vs1, Ts)

! FIXME:  this can cause problems if a worker gets data from all lower altitudes.  
if (maxval(Ts(i1:k1, i2:k2, i3:k3, i4:k4)) < 100) &
  call error_stop (dump_filename, 'input:plasma: too cold maximum Ts', ns, vs1, Ts)

end subroutine check_finite_plasma


subroutine check_finite_pertub(out_dir, t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

character(*), intent(in) :: out_dir
real(wp), intent(in) :: t_elapsed
integer, intent(in) :: worker_id

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn, vn1, vn2, vn3

character(:), allocatable :: dump_filename
character(8) :: wid

write(wid, '(I0)') worker_id
dump_filename = out_dir // "/dump_nonfinite_perturb_worker_" // trim(wid) // ".h5"

if (.not.all(ieee_is_finite(nn))) &
  call error_stop(dump_filename, 'perturb: non-finite Nn', t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

if (.not.all(ieee_is_finite(Tn))) &
  call error_stop(dump_filename, 'perturb: non-finite Tn', t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

if (.not.all(ieee_is_finite(vn1))) &
  call error_stop(dump_filename, 'perturb: non-finite vn1', t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

if (.not.all(ieee_is_finite(vn2))) &
  call error_stop(dump_filename, 'perturb: non-finite vn2', t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

if (.not.all(ieee_is_finite(vn3))) &
  call error_stop(dump_filename, 'perturb: non-finite vn3', t_elapsed, worker_id, nn, Tn, vn1, vn2, vn3)

end subroutine check_finite_pertub


subroutine check_finite_mag(out_dir, Br, Btheta, Bphi)

character(*), intent(in) :: out_dir
real(wp), dimension(:), intent(in) :: Br, Btheta, Bphi

character(:), allocatable :: dump_filename

dump_filename = out_dir // "/dump_nonfinite_mag.h5"

if (.not.all(ieee_is_finite(Br))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite Br', Br, Btheta, Bphi)
if (.not.all(ieee_is_finite(Btheta))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite Btheta', Br, Btheta, Bphi)
if (.not.all(ieee_is_finite(Bphi))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite Bphi', Br, Btheta, Bphi)

end subroutine check_finite_mag


subroutine check_finite_current(out_dir, J1,J2,J3)

character(*), intent(in) :: out_dir
real(wp), intent(in), dimension(:,:,:) :: J1, J2, J3

integer :: i1, k1, i2, k2, i3, k3

character(:), allocatable :: dump_filename

dump_filename = out_dir // "/dump_nonfinite_current.h5"


call ghost_bound(J1, i1,k1, i2,k2, i3,k3)

if (.not.all(ieee_is_finite(J1(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite J1', J1, J2, J3)

if (.not.all(ieee_is_finite(J2(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite J2', J1, J2, J3)

if (.not.all(ieee_is_finite(J3(i1:k1, i2:k2, i3:k3)))) &
  call error_stop(dump_filename, 'check_finite_mag: non-finite J3', J1, J2, J3)

end subroutine check_finite_current


end module sanity_check
