submodule (io) mag

use sanity_check, only : check_finite_mag
use pathlib, only : copyfile, mkdir

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
call mkdir(outdir//'/magfields/')
! call mkdir(outdir//'/magfields/input/')
! call copyfile(fieldpointfile, outdir//'/magfields/input/magfieldpoints.dat')

end procedure create_outdir_mag


module procedure output_magfields

select case (out_format)
case ('dat')
  call output_magfields_raw(outdir,ymd,UTsec,Br,Btheta,Bphi)
case ('h5')
  call output_magfields_hdf5(outdir,ymd,UTsec,Br,Btheta,Bphi)
case ('nc')
  call output_magfields_nc4(outdir,ymd,UTsec,Br,Btheta,Bphi)
case default
  error stop 'mag:output_magfields: unknown file format' // out_format
end select

call check_finite_mag(outdir, Br, Btheta, Bphi)


end procedure output_magfields


end submodule mag
