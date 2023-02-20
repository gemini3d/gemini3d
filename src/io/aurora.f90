submodule (io) io_aurora

use filesystem, only : mkdir, is_dir
use timeutils, only : date_filename
use mpimod, only : gather_send

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none (type, external)

interface ! aurora_*.f90

module subroutine output_aur_root_hdf5(filename,flagglow,iver)
character(*), intent(in) :: filename
integer, intent(in) :: flagglow
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine

end interface

contains


module procedure output_aur

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

character(:), allocatable :: outdir_aur

outdir_aur = outdir // '/aurmaps'

if(.not. is_dir(outdir_aur)) call mkdir(outdir_aur)

if (mpi_cfg%myid == 0) then
  call output_aur_root(date_filename(outdir_aur, ymd, UTsec), flagglow,iver, out_format)
else
  call output_aur_workers(iver)
end if

end procedure output_aur


subroutine output_aur_root(stem, flagglow, iver, out_format)

character(*), intent(in) :: stem, out_format
integer, intent(in) :: flagglow
real(wp), dimension(:,:,:), intent(in) :: iver

select case (out_format)
case ('h5')
  call output_aur_root_hdf5(stem // ".h5", flagglow, iver)
case default
  error stop 'ERROR:aurora:output_aur_root: unknown grid format' // out_format
end select

if(.not. all(ieee_is_finite(iver))) error stop 'ERROR: iverout: non-finite value(s)'

end subroutine output_aur_root


subroutine output_aur_workers(iver)
real(wp), dimension(:,:,:), intent(in) :: iver
!! SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)
!! The mpi'd dimensions are 2 and 3 so lwave needs to be permuted
!! to the first dimension for the canned routines to work.

integer :: iwave
real(wp), dimension(1:lx2,1:lx3) :: emistmp

do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_send(emistmp,tag%Aur)
end do


end subroutine output_aur_workers


end submodule io_aurora
