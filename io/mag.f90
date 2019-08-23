submodule (io) mag
implicit none
contains

module procedure create_outdir_mag
! subroutine create_outdir_mag(outdir,fieldpointfile)
!! CREATES OUTPUT DIRECTORY FOR MAGNETIC FIELD CALCULATIONS

integer :: ierr

!NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
ierr = mkdir(outdir//'/magfields/')
ierr = mkdir(outdir//'/magfields/input/')
ierr = copyfile(fieldpointfile, outdir//'/magfields/input/magfieldpoints.dat')

if (ierr/=0) error stop 'could not create magfields output directory'

end procedure create_outdir_mag


module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

character(:), allocatable :: outdir_composite, filenamefull
integer :: u


!FORM THE INPUT FILE NAME
outdir_composite=outdir//'/magfields/'
filenamefull=date_filename(outdir_composite,ymd,UTsec)
print *, '  Output file name (magnetic fields):  ',filenamefull
open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')

!> DUMP THE OUTPUT DATA
write(u) Br,Btheta,Bphi
close(u)

end procedure output_magfields

end submodule mag