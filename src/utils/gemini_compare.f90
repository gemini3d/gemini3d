program gemini_compare
!! for use from terminal/CMake
!! compares two directories that should have identical data
!! e.g. for CI

use compare_h5, only : check_plasma_hdf5
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none (type, external)

integer :: i, lx1, lx2all, lx3all, argc
character(1000) :: buf
character(:), allocatable :: new_path, ref_path
logical :: exists, all_ok

argc = command_argument_count()
if(argc /= 2) error stop 'specify top-level simulation path and reference path'

call get_command_argument(1, buf)
new_path = trim(buf)

call get_command_argument(2, buf)
ref_path = trim(buf)

call check_plasma_hdf5(new_path, ref_path, all_ok)


if (all_ok) stop "OK: gemini3d.compare: " // new_path // " == " // ref_path

error stop "gemini3d.compare FAIL " // new_path // " != " // ref_path

end program
