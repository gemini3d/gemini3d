submodule (io:mag) mag_hdf5

use timeutils, only: date_filename
use h5fortran, only: hdf5_file
implicit none

contains

module procedure output_magfields_hdf5
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(hdf5_file) :: hout

character(:), allocatable :: filenamefull

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.h5'
print *, '  Output file name (magnetic fields):  ',filenamefull

call hout%initialize(filenamefull, status='unknown',action='rw', comp_lvl=1)

call hout%write('/magfields/Br',     real(Br))
call hout%write('/magfields/Btheta', real(Btheta))
call hout%write('/magfields/Bphi',   real(Bphi))

call hout%finalize()

if(.not. all(ieee_is_finite(Br))) error stop 'Br: non-finite value(s)'
if(.not. all(ieee_is_finite(Btheta))) error stop 'Btheta: non-finite value(s)'
if(.not. all(ieee_is_finite(Bphi))) error stop 'Bphi: non-finite value(s)'

end procedure output_magfields_hdf5

end submodule mag_hdf5
