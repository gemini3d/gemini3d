submodule(io:io_aurora) io_aurora_hdf5

use, intrinsic:: iso_fortran_env, only: real32
use timeutils, only : date_filename
use hdf5_interface, only: hdf5_file

contains

module procedure output_aur_root
! subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver,zxden)
!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

type(hdf5_file) :: h5f

real(wp), dimension(1:lwave,1:lx2,1:lx3) :: ivertmp
real(wp), dimension(1:lwave,1:lx2all,1:lx3all) :: iverall

real(wp), dimension(1:lx2, 1:lx3) :: tmp3d
!! single emission subgrid

real(wp), dimension(1:lx2all, 1:lx3all) :: emisall
!real(real32) :: zxdenall(1:lx2all,1:lx3all,12,1:lx1)    !MZ - almost all of my arrays have x1 as the first dimension, suggest (lx1,lx2all,lx3all) to avoid heartburn later...
real(real32) :: zxdenall(1:lx1, 1:lx2all, 1:lx3all, 1:12)    !MZ - different dimension ordering so mpi works like we expect...
!! single emission total grid

real(wp), dimension(1:lx2all,1:lx3all,1:lwave) :: iverout  !< output array in the order scripts expect
integer :: i

character(:), allocatable :: outdir_composite, filenamefull, fstatus
integer :: u
logical :: exists

!! gather output from workers
do i=1,lwave
  tmp3d = iver(:,:,i)    !MZ - semantically confusing as this is a 2D array, but syntactically correct I think...
  call gather_recv(tmp3d, tagAur, emisall)
  iverout(:,:,i)=emisall
end do

excitationMPI : block

!real(real32) :: tmp4d(1:lx2, 1:lx3, 12, 1:lx1)  ! not used, just for notes
!do i=1,12
!  tmp4d=zxden(:,:,i,:)    !MZ - this is 4D, shouldn't it be 3D???
!  call gather_recv(tmp4d, tagZxden, zxdenall)    !MZ - this is going to call gather_recv4D_23, which assumes ghost cells (there are none) and also assumes the last dimension is species (looped over) and the second two dimensions are x2,x3 (mpi'd).  Also this needs to go into a 3D receive buffer, which then gets assigned to zxdenall(:,:,:,iwave)
!end do

  real(real32) :: tmpzxper(1:lx1, 1:lx2, 1:lx3)
  real(real32) :: tmpzxall(1:lx1, 1:lx2all, 1:lx3all)
  real(real32) :: tmpzx(1:lx2, 1:lx3, 1:lx1)

  do i=1,12  !< one species at a time
    tmpzx = zxden(:,:,i,:)
    !MZ - need lx1 to be first dim for mpi, maybe this needs to be done at the data source when GLOW get called???
    tmpzxper = reshape(zxden, [lx1,lx2,lx3], order=[3,1,2])

    call gather_recv(tmpzxper, tagZxden, tmpzxall)  !< 3D in, 3D+ghost out
    zxdenall(:,:,:,i)=tmpzxall     !assign buffer data to permanent 4D array
  end do
end block excitationMPI

!! create an output file
outdir_composite=outdir//'/aurmaps/'
filenamefull=date_filename(outdir_composite,ymd,UTsec) // '.h5'
!inquire(file=filenamefull, exist=exists)
!if (exists) then
!  fstatus = 'unknown'
!else
  print *, 'Output file name (auroral maps):  ',filenamefull
  fstatus = 'new'
!endif
call h5f%initialize(filenamefull, status=fstatus,action='rw',comp_lvl=1)

!! write data to file
if(flagswap/=1) then
  call h5f%add('/aurora/iverout', iverout)
  call h5f%add('/aurora/zxden', zxdenall)
else
  call h5f%add('/aurora/iverout', reshape(iverout,[lx3all,lx2all,lwave], order=[2,1,3]))
  call h5f%add('/aurora/zxden', reshape(zxdenall,[lx3all,lx2all,12,lx1], order=[2,1,3,4]))
end if

call h5f%finalize()

end procedure output_aur_root

end submodule io_aurora_hdf5
