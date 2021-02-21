program test_msis

use, intrinsic :: iso_fortran_env, only : real32

use h5fortran, only : hdf5_file
use assert, only : assert_allclose

implicit none (type, external)

integer :: argc
character(1000) :: argv
character(:), allocatable :: fnew, fref
real(real32) :: Dnew(9), Dref(9), Tnew(2), Tref(2), altnew, altref

argc = command_argument_count()
if(argc /= 2) error stop "./test_msis new_file ref_file"

call get_command_argument(1, argv)
fnew = trim(argv)

call get_command_argument(2, argv)
fref = trim(argv)

call reader(fnew, altnew, Dnew, Tnew)
call reader(fref, altref, Dref, Tref)

call assert_allclose(Tnew, Tref, rtol=1e-5, err_msg='mismatch: Tn')

call assert_allclose(Dnew, Dref, rtol=1e-5, err_msg='mismatch: Dn')

contains


subroutine reader(file, alt, Dn, Tn)

character(*), intent(in) :: file
real(real32), intent(out) :: alt, Dn(9), Tn(2)

real(real32) :: buf(1,1,1)  !< a priori for test file
type(hdf5_file) :: hf

call hf%initialize(file, status="old", action="read")

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


call hf%finalize()

end subroutine reader


end program
