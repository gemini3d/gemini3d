submodule (io) mag_hdf5

use timeutils, only: date_filename
use h5fortran, only: hdf5_file
implicit none

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(hdf5_file) :: hout

character(:), allocatable :: filenamefull

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call hout%initialize(filenamefull, status='unknown',action='rw')

call hout%write('/magfields/Br', Br)
call hout%write('/magfields/Btheta', Btheta)
call hout%write('/magfields/Bphi', Bphi)

call hout%finalize()

end procedure output_magfields

end submodule mag_hdf5
