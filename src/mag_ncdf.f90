submodule (io) mag_hdf5

use timeutils, only: date_filename
use nc4fortran, only: netcdf_file
implicit none

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(netcdf_file) :: hout

character(:), allocatable :: filenamefull
integer :: ierr

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call hout%initialize(filenamefull, ierr, status='unknown',action='rw',comp_lvl=1)

call hout%write('/magfields/Br', Br, ierr)
call hout%write('/magfields/Btheta', Btheta, ierr)
call hout%write('/magfields/Bphi', Bphi, ierr)

call hout%finalize(ierr)

end procedure output_magfields

end submodule mag_hdf5
