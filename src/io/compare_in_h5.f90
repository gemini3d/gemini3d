submodule (compare_h5) compare_in_h5

implicit none (type, external)

contains

module procedure check_plasma_input_hdf5

integer :: lx1, lx2all, lx3all, bad, i
character(:), allocatable :: new_file, ref_file

type(gemini_cfg) :: ref_cfg, cfg

character(6), parameter :: vars(3) = [character(6) :: "nsall", "Tsall", "vs1all"]

!> check sizes
call check_simsize(new_path, ref_path, lx1, lx2all, lx3all)

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
do i = 1,size(vars)

  bad = check_var(vars(i), new_file, ref_file, rtol, atol, lx1, lx2all, lx3all)

end do

all_ok = bad == 0


end procedure check_plasma_input_hdf5


integer function check_var(name, new_file, ref_file, rtol, atol, lx1, lx2all, lx3all) result(bad)

character(*), intent(in) :: new_file, ref_file
character(*), intent(in) :: name
real(wp), intent(in) :: rtol, atol
integer, intent(in) :: lx1, lx2all, lx3all
integer, parameter :: lsp = 7

type(hdf5_file) :: hnew, href

real, dimension(:,:,:,:), allocatable :: new, ref

bad = 0

allocate(new(lx1, lx2all, lx3all, lsp), ref(lx1, lx2all, lx3all, lsp))

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

call hnew%read(name, new)
call href%read(name, ref)

call hnew%finalize()
call href%finalize()

if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // name
if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(new_file) // " " // name

if(all(isclose(ref, new, rtol, atol))) then
  print '(A)', "OK: input: " // name
  return
endif

bad = 1

write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // name, &
' max diff:', maxval(abs(ref - new)), ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(ref)


end function check_var

end submodule compare_in_h5
