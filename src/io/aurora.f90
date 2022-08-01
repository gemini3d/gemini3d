submodule (io) io_aurora

use filesystem, only : mkdir
use h5fortran, only : hdf5_file
use timeutils, only : date_filename

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none (type, external)

contains


module procedure output_aur

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)

character(:), allocatable :: outdir_aur
type(hdf5_file) :: h
character(:), allocatable :: filename
integer :: ix2(2), ix3(2)

outdir_aur = outdir // '/aurmaps'

call mkdir(outdir_aur)

!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

!! output array in the order scripts expect

filename = date_filename(outdir_aur, ymd, UTsec) // ".h5"

print *, 'GEMINI3D:write aurora:  ',filename

call h%open(filename, action='w', comp_lvl=comp_lvl, mpi=.true.)

ix2 = [mpi_cfg%myid*lx2+1, (mpi_cfg%myid+1)*lx2]
ix3 = [mpi_cfg%myid*lx3+1, (mpi_cfg%myid+1)*lx3]
if(lx2 == 1) then
  ix2 = [1, 1]
elseif(lx3 == 1) then
  ix3 = [1, 1]
endif

! print *, "TRACE:gemini3d:output_aur: lx2all,lx3all,lwave,lx2,lx3",lx2all,lx3all,lwave,lx2,lx3
! print *, "TRACE:gemini3d:output_aur: myid: ", mpi_cfg%myid, " istart: ", [ix2(1), ix3(1), 1], " iend: ", [ix2(2), ix3(2), lwave]

call h%write('/aurora/iverout', real(iver), &
  dset_dims=[lx2all, lx3all, lwave], &
  istart=[ix2(1), ix3(1), 1], &
  iend=[ix2(2), ix3(2), lwave])

call h%close()

if(.not. all(ieee_is_finite(iver))) error stop 'ERROR:GEMINI3D:output_aur: non-finite value(s)'

end procedure output_aur


end submodule io_aurora
