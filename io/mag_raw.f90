submodule (io) mag_raw

use timeutils, only: date_filename

contains

module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

character(:), allocatable :: outdir_composite, filenamefull


!FORM THE INPUT FILE NAME
outdir_composite=outdir//'/magfields/'
filenamefull=date_filename(outdir_composite,ymd,UTsec) // '.dat'
print *, '  Output file name (magnetic fields):  ',filenamefull
block
  integer :: u
  open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')

  !> DUMP THE OUTPUT DATA
  write(u) Br,Btheta,Bphi
  close(u)
end block

if(.not. all(ieee_is_finite(Br))) error stop 'Br: non-finite value(s)'
if(.not. all(ieee_is_finite(Btheta))) error stop 'Btheta: non-finite value(s)'
if(.not. all(ieee_is_finite(Bphi))) error stop 'Bphi: non-finite value(s)'

end procedure output_magfields

end submodule mag_raw