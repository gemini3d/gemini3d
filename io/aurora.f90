submodule (io) io_aurora
implicit none
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
!subroutine output_aur(outdir,flagglow,ymd,UTsec,iver)

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

if (myid/=0) then
  call output_aur_workers(iver)
else
  call output_aur_root(outdir,flagglow,ymd,UTsec,iver)
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
  call gather_send(emistmp,tagAur)
end do


end procedure output_aur_workers


module procedure output_aur_root
! subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver)
!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

!real(wp), dimension(1:lx2,1:lwave,1:lx3) :: ivertmp
!real(wp), dimension(1:lx2all,1:lwave,1:lx3all) :: iverall
real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
real(wp), dimension(1:lwave,1:lx2all,1:lx3all) :: iverall

real(wp), dimension(1:lx2,1:lx3) :: emistmp                !single emission subgrid
real(wp), dimension(1:lx2all,1:lx3all) :: emisall          !single emission total grid
real(wp), dimension(1:lx2all,1:lx3all,1:lwave) :: iverout  !output array in the order scripts expect


character(:), allocatable :: outdir_composite, filenamefull
integer :: u
integer :: iwave

!!ivertmp=reshape(iver,[lx2,lwave,lx3],order=[1,3,2])
!ivertmp=reshape(iver,[lwave,lx2,lx3],order=[3,1,2])
!call gather_recv(ivertmp,tagAur,iverall)
do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_recv(emistmp,tagAur,emisall)
  iverout(:,:,iwave)=emisall
end do


!FORM THE INPUT FILE NAME
outdir_composite=outdir//'/aurmaps/'

filenamefull=date_filename(outdir_composite,ymd,UTsec)

print *, '  Output file name (auroral maps):  ',filenamefull
open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')

if(flagswap/=1) then
!  write(u) reshape(iverall,[lx2all,lwave,lx3all],order=[2,1,3])
!  write(u) reshape(iverall,[lx2all,lx3all,lwave],order=[2,3,1])
  write(u) iverout
else
!  write(u) reshape(iverall,[lx3all,lx2all,lwave],order=[3,2,1])
  write(u) reshape(iverout,[lx3all,lx2all,lwave],order=[2,1,3])
end if

close(u)
end procedure output_aur_root

end submodule io_aurora
