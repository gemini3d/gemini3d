submodule (io) mag

use sanity_check, only : check_finite_mag

implicit none (type, external)

interface ! mag_*.f90

module subroutine output_magfields_hdf5(outdir,ymd,UTsec,Br,Btheta,Bphi)
character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:), intent(in)  :: Br,Btheta,Bphi
end subroutine output_magfields_hdf5

module subroutine output_magfields_nc4(outdir,ymd,UTsec,Br,Btheta,Bphi)
character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:), intent(in)  :: Br,Btheta,Bphi
end subroutine output_magfields_nc4

module subroutine output_magfields_raw(outdir,ymd,UTsec,Br,Btheta,Bphi)
character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:), intent(in)  :: Br,Btheta,Bphi
end subroutine output_magfields_raw


end interface

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


module procedure output_magfields

select case (out_format)
case ('raw')
  call output_magfields_raw(outdir,ymd,UTsec,Br,Btheta,Bphi)
case ('h5', 'hdf5')
  call output_magfields_hdf5(outdir,ymd,UTsec,Br,Btheta,Bphi)
case ('nc', 'nc4')
  call output_magfields_nc4(outdir,ymd,UTsec,Br,Btheta,Bphi)
case default
  write(stderr,*) 'mag:output_magfields: unknown file format' // out_format
  error stop 6
end select

call check_finite_mag(Br, Btheta, Bphi)


end procedure output_magfields


end submodule mag
