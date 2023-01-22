submodule (io) io_cond
!! output conductivity to file
use filesystem, only : is_dir, mkdir
use timeutils, only : date_filename
use mpimod, only : gather_send

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none (type, external)

interface !< cond_*.f90
module subroutine output_cond_root_hdf5(filename, sig0, sigP, sigH)
character(*), intent(in) :: filename
real(wp), dimension(:,:,:), intent(in) :: sig0, sigP, sigH
end subroutine output_cond_root_hdf5
end interface

contains


module procedure output_cond

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

character(:), allocatable :: outdir_cond

outdir_cond = outdir // '/conductivity'

if(.not. is_dir(outdir_cond)) call mkdir(outdir_cond)

if (mpi_cfg%myid == 0) then
  call output_cond_root(date_filename(outdir, ymd, UTsec), sig0, sigP, sigH, out_format)
else
  call output_cond_workers(sig0, sigP, sigH)
end if

end procedure output_cond


subroutine output_cond_root(stem, sig0, sigP, sigH, out_format)
character(*), intent(in) :: stem, out_format
real(wp), dimension(:,:,:), intent(in) :: sig0, sigP, sigH

select case (out_format)
case ('h5')
  call output_cond_root_hdf5(stem // ".h5", sig0, sigP, sigH)
case default
  error stop 'io:cond:output_cond_root: unknown grid format' // out_format
end select

if(.not. all(ieee_is_finite(sig0))) error stop 'io:output_cond: non-finite sig0'
if(.not. all(ieee_is_finite(sigP))) error stop 'io:output_cond: non-finite sigP'
if(.not. all(ieee_is_finite(sigH))) error stop 'io:output_cond: non-finite sigH'

end subroutine output_cond_root


subroutine output_cond_workers(sig0, sigP, sigH)
!! SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.

real(wp), dimension(:,:,:), intent(in) :: sig0, sigP, sigH

call gather_send(sig0, tag%io_sig0)
call gather_send(sigP, tag%io_sigP)
call gather_send(sigH, tag%io_sigH)

end subroutine output_cond_workers


end submodule io_cond
