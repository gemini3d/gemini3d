program test_simsixe_nc4

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use nc4fortran, only: netcdf_file

implicit none (type, external)

logical :: exists

character(256) :: buf
character(:), allocatable :: path

type(netcdf_file) :: hf

integer :: lx(3)

call get_command_argument(1, buf)

path = trim(buf)

call hf%open(path, action='r')

call hf%read("lx", lx)

print *,lx

call hf%close()

end program
