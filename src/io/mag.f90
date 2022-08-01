submodule (io) mag

use sanity_check, only : check_finite_mag
use filesystem, only : mkdir
use timeutils, only : date_filename
use h5fortran, only : hdf5_file

implicit none (type, external)

contains

module procedure create_outdir_mag
! subroutine create_outdir_mag(outdir,fieldpointfile)
!! CREATES OUTPUT DIRECTORY FOR MAGNETIC FIELD CALCULATIONS

integer :: ierr

!NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
call mkdir(outdir // '/magfields')

end procedure create_outdir_mag


module procedure output_magfields
!! WE ASSUME THE ROOT PROCESS HAS ALREADY REDUCED THE MAGNETIC FIELD DATA

character(:), allocatable :: filename
type(hdf5_file) :: h

filename = date_filename(outdir // '/magfields', ymd, UTsec) // '.h5'
if(mpi_cfg%myid == 0) print *, '  Output file name (magnetic fields):  ',filename

call h%open(filename, action='rw', comp_lvl=comp_lvl, mpi=.false.)

call h%write('/magfields/Br',     real(Br))
call h%write('/magfields/Btheta', real(Btheta))
call h%write('/magfields/Bphi',   real(Bphi))

call h%close()

call check_finite_mag(outdir, Br, Btheta, Bphi)

end procedure output_magfields


end submodule mag
