program gemini_compare
!! for use from terminal/CMake
!! compares two directories that should have identical data
!! e.g. for CI

use compare_h5, only : check_plasma_output_hdf5, check_plasma_input_hdf5, check_grid
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none (type, external)

integer :: i,j, lx1, lx2all, lx3all, argc
character(1000) :: buf
character(6) :: which
character(:), allocatable :: new_path, ref_path
logical :: exists, all_ok
character(*), parameter :: help = './gemini3d.compare new_dir ref_dir [-which in|out]'

argc = command_argument_count()
if(argc < 2) error stop help

call get_command_argument(1, buf, status=i)
if (i/=0) error stop help
new_path = trim(buf)

call get_command_argument(2, buf, status=i)
if (i/=0) error stop help
ref_path = trim(buf)

buf = ""
if (argc > 2) then
  do j = 3,argc
    call get_command_argument(j, buf, status=i)
    if(i/=0) error stop help

    select case (buf)
    case ('-which')
      call get_command_argument(j+1, which, status=i)
      if (i/=0) error stop help
    end select
  end do
endif


all_ok = .false.
i = 0
!! in case .not. which /in {in,out}

i = index(buf, "grid")
!! we always check grid, but just to not falsely trip no-test error
call check_grid(new_path, ref_path)


i = index(buf, "out")
if (len_trim(buf) == 0 .or. i > 0) then
  call check_plasma_output_hdf5(new_path, ref_path, all_ok)

  if (all_ok) then
    print '(A)', "OK: gemini3d.compare: output: " // new_path // " == " // ref_path
  else
    error stop "gemini3d.compare output: FAIL " // new_path // " != " // ref_path
  endif
endif

i = index(buf, "in")
if (len_trim(buf) == 0 .or. i > 0) then
  call check_plasma_input_hdf5(new_path, ref_path, all_ok)

  if (all_ok) then
    print '(A)', "OK: gemini3d.compare: input: " // new_path // " == " // ref_path
  else
    error stop "gemini3d.compare input: FAIL " // new_path // " != " // ref_path
  endif
endif


if(len_trim(buf) > 0 .and. i < 1) error stop help
!! ensure at least one test run

print '(A)', "OK: gemini3d.compare: " // new_path

end program
