submodule(io:io_cond) io_cond_hdf5

use h5fortran, only: hdf5_file
use mpimod, only : gather_recv

implicit none (type, external)

contains

module procedure output_cond_root_hdf5
!! COLLECT COMPLETE DATA FROM WORKERS AND PROCESS FOR OUTPUT.
!! NO GHOST CELLS (I HOPE)

type(hdf5_file) :: hout

real(wp) :: tmp(1:lx1, 1:lx2, 1:lx3), tmpall(1:lx1, 1:lx2all, 1:lx3all)

!>  write file
print *, 'Output file name (conductivity):  ',filename

call hout%open(filename, action='w', comp_lvl=comp_lvl)

call gather_recv(tmp, tag%io_sig0, tmpall)
call hout%write('/sig0', tmpall)

call gather_recv(tmp, tag%io_sigP, tmpall)
call hout%write('/sigP', tmpall)

call gather_recv(tmp, tag%io_sigH, tmpall)
call hout%write('/sigH', tmpall)

call hout%close()

end procedure output_cond_root_hdf5

end submodule io_cond_hdf5
