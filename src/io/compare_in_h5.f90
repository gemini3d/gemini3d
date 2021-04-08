submodule (compare_h5) compare_in_h5

use pathlib, only : parent

implicit none (type, external)

contains

module procedure check_plasma_input_hdf5

logical :: dbug
integer :: bad
character(:), allocatable :: new_file, ref_file

type(gemini_cfg) :: ref_cfg, cfg

character(6), parameter :: var4(3) = [character(6) :: "nsall", "Tsall", "vs1all"]

dbug = .false.
if(present(debug)) dbug = debug

!> get input filename
ref_cfg%infile = new_path // '/inputs/config.nml'
ref_cfg%outdir = '.'  !< not used, just to pass checks
call read_configfile(ref_cfg)

cfg%infile = new_path // '/inputs/config.nml'
cfg%outdir = '.'  !< not used, just to pass checks
call read_configfile(cfg)

ref_file = ref_path // "/" // ref_cfg%indatfile
new_file = new_path // "/" // cfg%indatfile

!> check time
call check_time(new_file, ref_file)

!> check data
bad = 0

bad = bad + check_4d(new_file, ref_file, new_path, ref_path, var4, dbug)

check_plasma_input_hdf5 = bad == 0

end procedure check_plasma_input_hdf5


integer function check_4d(new_file, ref_file, new_path, ref_path, var4, debug)

character(*), intent(in) :: new_file, ref_file, new_path, ref_path
character(*), intent(in) :: var4(:)
logical, intent(in) :: debug

type(hdf5_file) :: href, hnew
integer :: i, lx1, lx2all, lx3all
real, allocatable :: new4(:,:,:,:), ref4(:,:,:,:)

call check_simsize(new_path, ref_path, lx1, lx2all, lx3all)

check_4d = 0

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

do i = 1,size(var4)

  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp))
  call hnew%read(var4(i), new4)
  call href%read(var4(i), ref4)

  if (.not.all(ieee_is_finite(ref4))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var4(i)
  if (.not.all(ieee_is_finite(new4))) error stop "NON-FINITE: " // file_name(new_file) // " " // var4(i)

  if(all(isclose(ref4, new4, rtol, atol))) then
    if(debug) print '(A)', "OK: input: " // var4(i)
  else
    check_4d = check_4d + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // var4(i), &
    ' max diff:', maxval(abs(ref4 - new4)), &
    ' max & min ref:', maxval(ref4), minval(ref4), ' max & min new:', maxval(new4), minval(new4)
  endif

  deallocate(new4, ref4)

end do

call hnew%finalize()
call href%finalize()

end function check_4d


end submodule compare_in_h5
