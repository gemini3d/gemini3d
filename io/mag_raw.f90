submodule (io) mag_raw

use timeutils, only: date_filename

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

character(:), allocatable :: outdir_composite, filenamefull
integer :: u


!FORM THE INPUT FILE NAME
outdir_composite=outdir//'/magfields/'
filenamefull=date_filename(outdir_composite,ymd,UTsec) // '.dat'
print *, '  Output file name (magnetic fields):  ',filenamefull
open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')

!> DUMP THE OUTPUT DATA
write(u) Br,Btheta,Bphi
close(u)

end procedure output_magfields

end submodule mag_raw