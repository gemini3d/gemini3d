submodule (io:mag) mag_hdf5

use timeutils, only: date_filename
use h5fortran, only: hdf5_file

implicit none (type, external)

contains

module procedure output_magfields_hdf5
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(hdf5_file) :: hout

character(:), allocatable :: filenamefull

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call hout%initialize(filenamefull, status='unknown',action='rw', comp_lvl=comp_lvl)

call hout%write('/magfields/Br',     real(Br))
call hout%write('/magfields/Btheta', real(Btheta))
call hout%write('/magfields/Bphi',   real(Bphi))

call hout%finalize()

end procedure output_magfields_hdf5

end submodule mag_hdf5
