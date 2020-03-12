submodule(io:io_aurora) io_aurora_hdf5

use timeutils, only : date_filename
use nc4fortran, only: netcdf_file
implicit none

contains

module procedure output_aur_root
! subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver)
!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

type(netcdf_file) :: hout

real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
real(wp), dimension(1:lwave,1:lx2all,1:lx3all) :: iverall

real(wp), dimension(1:lx2,1:lx3) :: emistmp                !< single emission subgrid
real(wp), dimension(1:lx2all,1:lx3all) :: emisall          !< single emission total grid
real(wp), dimension(1:lx2all,1:lx3all,1:lwave) :: iverout  !< output array in the order scripts expect
integer :: iwave

character(:), allocatable :: outdir_composite, filenamefull, fstatus
integer :: ierr

!! gather output from workers
do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_recv(emistmp,tagAur,emisall)
  iverout(:,:,iwave)=emisall
end do

!! create an output file
outdir_composite=outdir//'/aurmaps/'
filenamefull=date_filename(outdir_composite,ymd,UTsec) // '.h5'
print *, 'Output file name (auroral maps):  ',filenamefull
fstatus = 'new'
call hout%initialize(filenamefull, ierr, status=fstatus,action='rw',comp_lvl=1)

!! write data to file
if(flagswap/=1) then
  call hout%write('/aurora/iverout', iverout, ierr)
else
  call hout%write('/aurora/iverout', reshape(iverout,[lx3all,lx2all,lwave],order=[2,1,3]), ierr)
end if

call hout%finalize(ierr)

end procedure output_aur_root

end submodule io_aurora_hdf5