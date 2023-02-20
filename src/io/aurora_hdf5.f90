submodule(io:io_aurora) io_aurora_hdf5

use h5fortran, only: hdf5_file
use mpimod, only : gather_recv

implicit none (type, external)

contains

module procedure output_aur_root_hdf5
  !! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
  !! NO GHOST CELLS (I HOPE)

  type(hdf5_file) :: hout

  real(wp), dimension(1:lx2,1:lx3) :: emistmp                !< single emission subgrid
  real(wp), dimension(1:lx2all,1:lx3all) :: emisall          !< single emission total grid
  real(wp), dimension(1:lx2all,1:lx3all,1:lwave) :: iverout  !< output array in the order scripts expect
  integer :: iwave

  !! gather output from workers
  do iwave=1,lwave
    emistmp=iver(:,:,iwave)
    call gather_recv(emistmp,tag%Aur,emisall)
    iverout(:,:,iwave)=emisall
  end do

  !! create an output file
  print *, 'write aurora:  ',filename

  call hout%open(filename, action='w', comp_lvl=comp_lvl)

  !! write data to file
  call hout%write('/aurora/iverout', real(iverout))

  call hout%close()
end procedure output_aur_root_hdf5

end submodule io_aurora_hdf5
