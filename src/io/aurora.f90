submodule (io) io_aurora

use pathlib, only : mkdir

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none (type, external)

interface ! aurora_*.f90

module subroutine output_aur_root_raw(outdir,flagglow,ymd,UTsec,iver)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_root_raw

module subroutine output_aur_root_hdf5(outdir,flagglow,ymd,UTsec,iver)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_root_hdf5

module subroutine output_aur_root_nc4(outdir,flagglow,ymd,UTsec,iver)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_root_nc4

end interface

contains


subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver, out_format)
character(*), intent(in) :: outdir, out_format
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver

select case (out_format)
case ('dat')
  call output_aur_root_raw(outdir,flagglow,ymd,UTsec,iver)
case ('h5')
  call output_aur_root_hdf5(outdir,flagglow,ymd,UTsec,iver)
case ('nc')
  call output_aur_root_nc4(outdir,flagglow,ymd,UTsec,iver)
case default
  write(stderr,*) 'aurora:output_aur_root: unknown grid format' // out_format
  error stop 2
end select

if(.not. all(ieee_is_finite(iver))) error stop 'iverout: non-finite value(s)'

end subroutine output_aur_root

module procedure create_outdir_aur
!subroutine create_outdir_aur(outdir)
!! CREATES OUTPUT DIRECTORY FOR Auroral CALCULATIONS
integer :: ierr

!NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
ierr = mkdir(outdir//'/aurmaps/')

end procedure create_outdir_aur


module procedure output_aur
!subroutine output_aur(outdir,flagglow,ymd,UTsec,iver, out_format)

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

if (mpi_cfg%myid/=0) then
  call output_aur_workers(iver)
else
  call output_aur_root(outdir,flagglow,ymd,UTsec,iver, out_format)
end if

end procedure output_aur


module procedure output_aur_workers
!subroutine output_aur_workers(iver)
!! SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)
!! The mpi'd dimensions are 2 and 3 so lwave needs to be permuted
!! to the first dimension for the canned routines to work.

!real(wp), dimension(1:lx2,1:lwave,1:lx3) :: ivertmp
real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
integer :: iwave
real(wp), dimension(1:lx2,1:lx3) :: emistmp


!!ivertmp=reshape(iver,[lx2,lwave,lx3],order=[1,3,2])
!ivertmp=reshape(iver,[lwave,lx2,lx3],order=[3,1,2])
!
!!------- SEND AURORA PARAMETERS TO ROOT
!call gather_send(ivertmp,tagAur)

do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_send(emistmp,tag%Aur)
end do


end procedure output_aur_workers


end submodule io_aurora
