submodule (compare_h5) compare_in_h5

use pathlib, only : parent

implicit none (type, external)

contains

module procedure check_plasma_input_hdf5

logical :: dbug
integer :: bad
character(:), allocatable :: new_file, ref_file

type(gemini_cfg) :: ref_cfg, cfg

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

bad = bad + check_initcond(new_file, ref_file, new_path, ref_path, dbug)

bad = bad + check_precip(new_path, ref_path, cfg, dbug)

bad = bad + check_Efield(new_path, ref_path, cfg, dbug)

check_plasma_input_hdf5 = bad == 0

end procedure check_plasma_input_hdf5


integer function check_initcond(new_file, ref_file, new_path, ref_path, debug) result(bad)

character(*), intent(in) :: new_path, ref_path, new_file, ref_file
logical, intent(in) :: debug

character(6), parameter :: var(3) = [character(6) :: "nsall", "Tsall", "vs1all"]

type(hdf5_file) :: href, hnew
integer :: i, lx1, lx2all, lx3all

real, allocatable :: new4(:,:,:,:), ref4(:,:,:,:)

call check_simsize(new_path, ref_path, lx1, lx2all, lx3all)

bad = 0

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

do i = 1,size(var)


  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp))
  call hnew%read(var(i), new4)
  call href%read(var(i), ref4)

  if (.not.all(ieee_is_finite(ref4))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var(i)
  if (.not.all(ieee_is_finite(new4))) error stop "NON-FINITE: " // file_name(new_file) // " " // var(i)

  if(all(isclose(ref4, new4, rtol, atol))) then
    if(debug) print '(A)', "OK: input: " // var(i)
  else
    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:init_cond: " // file_name(new_file) // " " // var(i), &
    ' max diff:', maxval(abs(ref4 - new4)), &
    ' max & min ref:', maxval(ref4), minval(ref4), ' max & min new:', maxval(new4), minval(new4)
  endif

  deallocate(new4, ref4)

end do

call hnew%finalize()
call href%finalize()

end function check_initcond


integer function check_precip(new_path, ref_path, cfg, debug) result(bad)

character(*), intent(in) :: new_path, ref_path
class(gemini_cfg), intent(in) :: cfg
logical, intent(in) :: debug

character(3), parameter :: var(*) = [character(3) :: "Qp", "E0p"]

character(:),allocatable :: new_file, ref_file, suffix
type(hdf5_file) :: href, hnew
integer :: i, lx2, lx3, ymd(3)
real(wp) :: UTsec, t

real, allocatable :: new(:,:), ref(:,:)

suffix = get_suffix(cfg%indatsize)

call check_simsize2(new_path // "/" // cfg%precdir, ref_path // "/" // cfg%precdir, lx2, lx3)

bad = 0
t = 0
ymd = cfg%ymd0
UTsec = cfg%UTsec0

do while (t <= cfg%tdur)

  new_file = date_filename(new_path // "/" // cfg%precdir, ymd, UTsec) // suffix
  ref_file = date_filename(ref_path // "/" // cfg%precdir, ymd, UTsec) // suffix

  call hnew%initialize(new_file, status='old',action='r')
  call href%initialize(ref_file, status='old',action='r')

  do i = 1,size(var)

    allocate(new(lx2, lx3), ref(lx2, lx3))
    call hnew%read(var(i), new)
    call href%read(var(i), ref)

    if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var(i)
    if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(new_file) // " " // var(i)

    if(all(isclose(ref, new, rtol, atol))) then
      if(debug) print '(A)', "OK: input:precip " // var(i) // " " // new_file
    else
      bad = bad + 1

      write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:precip: " // file_name(new_file) // " " // var(i), &
      ' max diff:', maxval(abs(ref - new)), &
      ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(new)
    endif

    deallocate(new, ref)

  end do

  call hnew%finalize()
  call href%finalize()

  call dateinc(cfg%dtprec, ymd, UTsec)

  t = t + cfg%dtprec
end do


end function check_precip


integer function check_Efield(new_path, ref_path, cfg, debug) result(bad)

character(*), intent(in) :: new_path, ref_path
class(gemini_cfg), intent(in) :: cfg
logical, intent(in) :: debug

character(8), parameter :: var(*) = [character(8) :: "Exit", "Eyit", "Vminx1it", "Vmaxx1it"]

character(:),allocatable :: new_file, ref_file, suffix
type(hdf5_file) :: href, hnew
integer :: i, lx2, lx3, ymd(3)
real(wp) :: UTsec, t

real, allocatable :: new(:,:), ref(:,:)

suffix = get_suffix(cfg%indatsize)

call check_simsize2(new_path // "/" // cfg%E0dir, ref_path // "/" // cfg%E0dir, lx2, lx3)

bad = 0
t = 0
ymd = cfg%ymd0
UTsec = cfg%UTsec0

do while (t <= cfg%tdur)

  new_file = date_filename(new_path // "/" // cfg%E0dir, ymd, UTsec) // suffix
  ref_file = date_filename(ref_path // "/" // cfg%E0dir, ymd, UTsec) // suffix

  call hnew%initialize(new_file, status='old',action='r')
  call href%initialize(ref_file, status='old',action='r')
  do i = 1,size(var)


    allocate(new(lx2, lx3), ref(lx2, lx3))
    call hnew%read(var(i), new)
    call href%read(var(i), ref)

    if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var(i)
    if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(new_file) // " " // var(i)

    if(all(isclose(ref, new, rtol, atol))) then
      if(debug) print '(A)', "OK: input:precip " // var(i) // " " // new_file
    else
      bad = bad + 1

      write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:precip: " // file_name(new_file) // " " // var(i), &
      ' max diff:', maxval(abs(ref - new)), &
      ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(new)
    endif

    deallocate(new, ref)

  end do

  call hnew%finalize()
  call href%finalize()

  call dateinc(cfg%dtE0, ymd, UTsec)

  t = t + cfg%dtE0

end do


end function check_Efield


end submodule compare_in_h5
