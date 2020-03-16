program test_simsixe_nc4

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use nc4fortran, only: netcdf_file

implicit none

logical :: exists

character(256) :: buf
character(:), allocatable :: path

type(netcdf_file) :: hf

integer :: lx(3)

call get_command_argument(1, buf)

path = buf

inquire(file=path, exist=exists)
if (.not.exists) then
   write(stderr,'(A,/,A)') 'ERROR: reader_nc4:get_simsize3: generate grid with script--grid not present: ', path
   error stop 77
endif

call hf%initialize(path, status='old', action='r')

call hf%read("lx", lx)

print *,lx

call hf%finalize()

end program