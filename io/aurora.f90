submodule (io) io_aurora

implicit none

interface ! aurora_*.f90

module subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver,zxden)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
real(real32), intent(in) :: zxden(:,:,:,:)
end subroutine output_aur_root

end interface

contains


module procedure create_outdir_aur
!subroutine create_outdir_aur(outdir)
!! CREATES OUTPUT DIRECTORY FOR Auroral CALCULATIONS
integer :: ierr

!NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
ierr = mkdir(outdir//'/aurmaps/')
if (ierr /= 0) error stop 'could not create auroral map output directory'

end procedure create_outdir_aur


module procedure output_aur
!subroutine output_aur(outdir,flagglow,ymd,UTsec,iver,zxden)

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

if (myid/=0) then
  call output_aur_workers(iver,zxden)
else
  call output_aur_root(outdir,flagglow,ymd,UTsec,iver,zxden)
end if

end procedure output_aur


module procedure output_aur_workers
!subroutine output_aur_workers(iver, zxden)
!! SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)
!! The mpi'd dimensions are 2 and 3 so lwave needs to be permuted
!! to the first dimension for the canned routines to work.

!real(wp), dimension(1:lx2,1:lwave,1:lx3) :: ivertmp
!real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
integer :: iwave
real(wp), dimension(1:lx2,1:lx3) :: emistmp


!!ivertmp=reshape(iver,[lx2,lwave,lx3],order=[1,3,2])
!ivertmp=reshape(iver,[lwave,lx2,lx3],order=[3,1,2])
!
!!------- SEND AURORA PARAMETERS TO ROOT
!call gather_send(ivertmp,tagAur)

do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_send(emistmp,tagAur)
end do


end procedure output_aur_workers


end submodule io_aurora
