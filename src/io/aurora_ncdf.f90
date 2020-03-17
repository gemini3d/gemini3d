submodule(io:io_aurora) io_aurora_nc4

use timeutils, only : date_filename
use nc4fortran, only: netcdf_file
implicit none

contains

module procedure output_aur_root_nc4
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

!! gather output from workers
do iwave=1,lwave
  emistmp=iver(:,:,iwave)
  call gather_recv(emistmp,tag%Aur,emisall)
  iverout(:,:,iwave)=emisall
end do

!! create an output file
outdir_composite=outdir//'/aurmaps/'
filenamefull=date_filename(outdir_composite,ymd,UTsec) // '.h5'
print *, 'Output file name (auroral maps):  ',filenamefull
fstatus = 'new'
call hout%initialize(filenamefull, status=fstatus,action='rw',comp_lvl=1)

!! write data to file
if(flagswap/=1) then
  call hout%write('/aurora/iverout', iverout)
else
  call hout%write('/aurora/iverout', reshape(iverout,[lx3all,lx2all,lwave],order=[2,1,3]))
end if

call hout%finalize()

end procedure output_aur_root_nc4

end submodule io_aurora_nc4