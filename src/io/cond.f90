submodule (io) io_cond
!! output conductivity to file
use h5fortran, only : hdf5_file
use filesystem, only : mkdir
use timeutils, only : date_filename

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

implicit none (type, external)

contains


module procedure output_cond

!! A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
!! BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
!! VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)
!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

character(:), allocatable :: outdir_cond, filename
type(hdf5_file) :: h

integer, dimension(3) :: dims, i0, i1
integer, dimension(2) :: ix1, ix2, ix3

outdir_cond = outdir // '/conductivity'

call mkdir(outdir_cond)

filename = date_filename(outdir_cond, ymd, UTsec) // ".h5"

print *, 'Output file name (conductivity):  ',filename

dims = [lx1, lx2all, lx3all]

ix1 = [1, lx1]
ix2 = [mpi_cfg%myid*lx2+1, (mpi_cfg%myid+1)*lx2]
ix3 = [mpi_cfg%myid*lx3+1, (mpi_cfg%myid+1)*lx3]
if(lx2 == 1) then
  ix2 = [1, 1]
elseif(lx3 == 1) then
  ix3 = [1, 1]
endif

i0 = [ix1(1), ix2(1), ix3(1)]
i1 = [ix1(2), ix2(2), ix3(2)]

call h%open(filename, action='w', comp_lvl=comp_lvl, mpi=.true.)

call h%write('/sig0', real(sig0), dset_dims=dims, istart=i0, iend=i1)
call h%write('/sigP', real(sigP), dset_dims=dims, istart=i0, iend=i1)
call h%write('/sigH', real(sigH), dset_dims=dims, istart=i0, iend=i1)

call h%close()

if(.not. all(ieee_is_finite(sig0))) error stop 'ERROR:GEMINI3D:io:output_cond: non-finite sig0'
if(.not. all(ieee_is_finite(sigP))) error stop 'ERROR:GEMINI3D:io:output_cond: non-finite sigP'
if(.not. all(ieee_is_finite(sigH))) error stop 'ERROR:GEMINI3D:io:output_cond: non-finite sigH'

end procedure output_cond


end submodule io_cond
