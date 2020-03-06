submodule (io) mag

implicit none

contains

module procedure create_outdir_mag
! subroutine create_outdir_mag(outdir,fieldpointfile)
!! CREATES OUTPUT DIRECTORY FOR MAGNETIC FIELD CALCULATIONS

integer :: ierr

!NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
if ( mkdir(outdir//'/magfields/') /= 0 ) error stop 'could not create magfields output directory'
if ( mkdir(outdir//'/magfields/input/') /= 0 ) error stop 'could not create magfields/input directory'
if ( copyfile(fieldpointfile, outdir//'/magfields/input/magfieldpoints.dat') /=0 ) error stop 'could not create magfields dat file'

end procedure create_outdir_mag


end submodule mag
