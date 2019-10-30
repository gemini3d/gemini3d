submodule (io) mag_hdf5

use timeutils, only: date_filename
use hdf5_interface, only: hdf5_file

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(hdf5_file) :: h5f

character(:), allocatable :: filenamefull

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call h5f%initialize(filenamefull, status='unknown',action='rw',comp_lvl=1)

call h5f%add('/magfields/Br', Br)
call h5f%add('/magfields/Btheta', Btheta)
call h5f%add('/magfields/Bphi', Bphi)

call h5f%finalize()

end procedure output_magfields

end submodule mag_hdf5
