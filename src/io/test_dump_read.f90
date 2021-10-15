program test_dump_read
!! quick checks that dump file was written as expected.
use iso_c_binding, only : c_null_char
use h5fortran, only : hdf5_file

implicit none (type, external)

type(hdf5_file) :: h

character(1000) :: buf
character(:), allocatable :: filename, mode
logical :: exists

if (command_argument_count() /= 2) error stop "please give: filename expected_mode"

call get_command_argument(1, buf)
mode = trim(buf)

call get_command_argument(2, buf)
filename = trim(buf)

inquire(file=filename, exist=exists)
if(.not.exists) error stop filename // " does not exist"

call h%open(filename, action="r")
call h%read("/info", buf)
call h%close()

if (buf /= mode) error stop "dump mode mismatch: expected " // mode // " but got " // trim(buf)

end program
