program gen_msis_test_in
!! Generate test_data/msis/msis_test_in.h5
!! assumes directory to write data into already exists.

use, intrinsic :: iso_fortran_env, only : real32
use h5fortran, only : hdf5_file
implicit none (type, external)

integer :: i, msis_version
character(1000) :: buf
character(:), allocatable :: filename

integer, parameter :: doy = 50
real(real32), parameter :: &
Ap(7) = 5.0, &
f107a = 109.5, &
f107 = 109.9, &
UTsec = 18000.0, &
alt(1,1,1) = 300.0, &
glat(1,1,1) = 66.0, &
glon(1,1,1) = 210.0

type(hdf5_file) :: f

if(command_argument_count() /= 2) error stop "please specify: generated_filename msis_version"

call get_command_argument(1, buf, status=i)
if(i/=0) error stop "please specify file name to generate MSIS test input data"
filename = trim(buf)

call get_command_argument(2, buf, status=i)
read(buf,'(i2)') msis_version

call f%open(filename, action="w")

call f%write("/msis_version", msis_version)
call f%write("/doy", doy)
call f%write("/Ap", Ap)
call f%write("/f107", f107)
call f%write("/f107a", f107a)
call f%write("/UTsec", UTsec)
call f%write("/alt", alt)
call f%write("/glat", glat)
call f%write("/glon", glon)

call f%close()

end program
