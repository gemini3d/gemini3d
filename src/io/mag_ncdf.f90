submodule (io:mag) mag_nc4

use timeutils, only: date_filename
use nc4fortran, only: netcdf_file

implicit none (type, external)

contains

module procedure output_magfields_nc4
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

type(netcdf_file) :: hout

character(:), allocatable :: filenamefull

filenamefull = date_filename(outdir // '/magfields/',ymd,UTsec) // '.nc'
print *, '  Output file name (magnetic fields):  ',filenamefull

call hout%open(filenamefull, action='rw',comp_lvl=comp_lvl)

call hout%write('magfields/Br', Br)
call hout%write('magfields/Btheta', Btheta)
call hout%write('magfields/Bphi', Bphi)

call hout%close()

end procedure output_magfields_nc4

end submodule mag_nc4
