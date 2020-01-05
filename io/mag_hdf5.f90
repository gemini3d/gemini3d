submodule (io) mag_hdf5

use timeutils, only: date_filename
use h5fortran, only: hdf5_file

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(hdf5_file) :: h5f

character(:), allocatable :: filenamefull
integer :: ierr

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call h5f%initialize(filenamefull, ierr, status='unknown',action='rw',comp_lvl=1)

call h5f%write('/magfields/Br', Br, ierr)
call h5f%write('/magfields/Btheta', Btheta, ierr)
call h5f%write('/magfields/Bphi', Bphi, ierr)

call h5f%finalize(ierr)

end procedure output_magfields

end submodule mag_hdf5
