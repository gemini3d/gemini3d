submodule (io:mag) mag_raw

use timeutils, only: date_filename

implicit none (type, external)

contains

module procedure output_magfields_raw
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

end procedure output_magfields_raw

end submodule mag_raw