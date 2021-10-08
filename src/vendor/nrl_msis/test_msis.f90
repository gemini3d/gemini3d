program test_msis

use, intrinsic :: iso_fortran_env, only : real32

use h5fortran, only : hdf5_file
use assert, only : assert_isclose

implicit none (type, external)

integer :: argc, msis_version_new, msis_version_ref
character(1000) :: argv
character(:), allocatable :: fnew, fref
real(real32) :: Dnew(9), Dref(9), Tnew(2), Tref(2), altnew, altref

argc = command_argument_count()
if(argc /= 2) error stop "./test_msis new_file ref_file"

call get_command_argument(1, argv)
fnew = trim(argv)

call get_command_argument(2, argv)
fref = trim(argv)

call reader(fnew, msis_version_new, altnew, Dnew, Tnew)
call reader(fref, msis_version_ref, altref, Dref, Tref)

if(msis_version_new /= msis_version_ref) error stop "msis_version differs"
call assert_isclose(Tnew, Tref, rtol=1e-5, err_msg='mismatch: Tn')
call assert_isclose(Dnew, Dref, rtol=1e-5, err_msg='mismatch: Dn')

contains


subroutine reader(file, msis_version, alt, Dn, Tn)

character(*), intent(in) :: file
integer, intent(out) :: msis_version
real(real32), intent(out) :: alt, Dn(9), Tn(2)

real(real32) :: buf(1,1,1)  !< a priori for test file
type(hdf5_file) :: hf

call hf%open(file, action="r")

call hf%read("/msis_version", msis_version)

call hf%read("/alt", buf)
alt = buf(1,1,1)

call hf%read("/nHe", buf)
Dn(1) = buf(1,1,1)
call hf%read("/nO", buf)
Dn(2) = buf(1,1,1)
call hf%read("/nN2", buf)
Dn(3) = buf(1,1,1)
call hf%read("/nO2", buf)
Dn(4) = buf(1,1,1)
call hf%read("/nAr", buf)
Dn(5) = buf(1,1,1)
call hf%read("/TotalMassDensity", buf)
Dn(6) = buf(1,1,1)
call hf%read("/nH", buf)
Dn(7) = buf(1,1,1)
call hf%read("/nN", buf)
Dn(8) = buf(1,1,1)
call hf%read("/nOana", buf)
Dn(9) = buf(1,1,1)

call hf%read("/Tn", buf)
Tn(2) = buf(1,1,1)
call hf%read("/Texo", buf)
Tn(1) = buf(1,1,1)


call hf%close()

end subroutine reader


end program
