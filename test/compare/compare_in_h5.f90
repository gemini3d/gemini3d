submodule (compare_h5) compare_in_h5

implicit none (type, external)

contains

module procedure check_plasma_input_hdf5

integer :: bad

type(hdf5_file) :: hnew, href
type(gemini_cfg) :: ref_cfg, cfg

!> get input filename
ref_cfg%infile = new_path // '/inputs/config.nml'
ref_cfg%outdir = ref_path  !< not used, just to pass checks
call read_configfile(ref_cfg)
if(.not. is_file(ref_cfg%indatfile)) error stop "compare: ref initcond file not found: " // ref_cfg%indatfile

cfg%infile = new_path // '/inputs/config.nml'
cfg%outdir = new_path  !< not used, just to pass checks
call read_configfile(cfg)
if(.not. is_file(cfg%indatfile)) error stop "compare: new initcond file not found: " // cfg%indatfile

! print '(a)', "check_plasma_input_hdf5: opening reference " // cfg%indatfile
call href%open(ref_cfg%indatfile, action='r')

! print '(a)', "check_plasma_input_hdf5: opening data " // cfg%indatfile
call hnew%open(cfg%indatfile, action='r')

!> check time
call check_time(hnew, href)

!> check data
bad = 0

bad = bad + check_initcond(hnew, href, new_path, ref_path, P)

call hnew%close()
call href%close()

if (cfg%flagprecfile == 1) then
  bad = bad + check_precip(ref_cfg, cfg, P)
endif

if (cfg%flagE0file == 1) then
  bad = bad + check_Efield(ref_cfg, cfg, P)
endif

check_plasma_input_hdf5 = bad == 0

end procedure check_plasma_input_hdf5


integer function check_initcond(hnew, href, new_path, ref_path, P) result(bad)

type(hdf5_file), intent(in) :: hnew, href
character(*), intent(in) :: new_path, ref_path
class(params), intent(in) :: P

character(6), parameter :: var(3) = [character(6) :: "nsall", "Tsall", "vs1all"]

integer :: i, lx1, lx2all, lx3all

real, allocatable :: new4(:,:,:,:), ref4(:,:,:,:)

call check_simsize(new_path, ref_path, lx1, lx2all, lx3all)

bad = 0


do i = 1,size(var)

  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp))
  call hnew%read(var(i), new4)
  call href%read(var(i), ref4)

  if (.not.all(ieee_is_finite(ref4))) error stop "NON-FINITE: " // file_name(href%filename) // " " // var(i)
  if (.not.all(ieee_is_finite(new4))) error stop "NON-FINITE: " // file_name(hnew%filename) // " " // var(i)

  if(all(isclose(ref4, new4, real(rtol), real(atol)))) then
    if(P%debug) print '(A)', "OK: input: " // var(i)
  else
    bad = bad + 1

    write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:init_cond: " // file_name(hnew%filename) // " " // var(i), &
    ' max diff:', maxval(abs(ref4 - new4)), &
    ' max & min ref:', maxval(ref4), minval(ref4), ' max & min new:', maxval(new4), minval(new4)
  endif

  deallocate(new4, ref4)

end do

if(bad /= 0) call plot_diff(hnew%filename, href%filename, "init_cond", "in", P)

print '(A)', "OK: input: " // hnew%filename

end function check_initcond


integer function check_precip(ref, new, P) result(bad)

class(gemini_cfg), intent(in) :: ref, new
class(params), intent(in) :: P

character(3), parameter :: var(*) = [character(3) :: "Qp", "E0p"]

character(:),allocatable :: new_file, ref_file
type(hdf5_file) :: href, hnew
integer :: i, lx2, lx3, ymd(3)
real(wp) :: UTsec, t

real, allocatable :: Anew(:,:), Aref(:,:)

call check_simsize2(new%precdir, ref%precdir, lx2, lx3)

bad = 0
t = 0
ymd = new%ymd0
UTsec = new%UTsec0

do while (t <= new%tdur)

  new_file = date_filename(new%precdir, ymd, UTsec) // ".h5"
  ref_file = date_filename(ref%precdir, ymd, UTsec) // ".h5"

  if(.not. is_file(new_file)) error stop "compare:check_precip: new precip file not found: " // new_file
  if(.not. is_file(ref_file)) error stop "compare:check_precip: ref precip file not found: " // ref_file

  call hnew%open(new_file, action='r')
  call href%open(ref_file, action='r')

  do i = 1,size(var)

    allocate(Anew(lx2, lx3), Aref(lx2, lx3))
    call hnew%read(var(i), Anew)
    call href%read(var(i), Aref)

    if (.not.all(ieee_is_finite(Aref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var(i)
    if (.not.all(ieee_is_finite(Anew))) error stop "NON-FINITE: " // file_name(new_file) // " " // var(i)

    if(all(isclose(Aref, Anew, real(rtol), real(atol)))) then
      if(P%debug) print '(A)', "OK: input:precip " // var(i) // " " // new_file
    else
      bad = bad + 1

      write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:precip: " // file_name(new_file) // " " // var(i), &
      ' max diff:', maxval(abs(Aref - Anew)), &
      ' max & min ref:', maxval(Aref), minval(Aref), ' max & min new:', maxval(Anew), minval(Anew)
    endif

    deallocate(Anew, Aref)

  end do

  call hnew%close()
  call href%close()

  call dateinc(new%dtprec, ymd, UTsec)

  t = t + new%dtprec
end do

if (bad == 0) print '(a,1x,a)', "OK: precip: ", new%precdir

end function check_precip


integer function check_Efield(ref, new, P) result(bad)

class(gemini_cfg), intent(in) :: ref, new
class(params), intent(in) :: P

character(8), parameter :: var(*) = [character(8) :: "Exit", "Eyit", "Vminx1it", "Vmaxx1it"]

character(:),allocatable :: new_file, ref_file
type(hdf5_file) :: href, hnew
integer :: i, lx2, lx3, ymd(3)
real(wp) :: UTsec, t

real, allocatable :: Anew(:,:), Aref(:,:)

call check_simsize2(new%E0dir, new%E0dir, lx2, lx3)

bad = 0
t = 0
ymd = new%ymd0
UTsec = new%UTsec0

do while (t <= new%tdur)

  new_file = date_filename(new%E0dir, ymd, UTsec) // ".h5"
  ref_file = date_filename(ref%E0dir, ymd, UTsec) // ".h5"

  if(.not. is_file(new_file)) error stop "compare:check_Efield: new Efield file not found: " // new_file
  if(.not. is_file(ref_file)) error stop "compare:check_Efield: ref Efield file not found: " // ref_file

  call hnew%open(new_file, action='r')
  call href%open(ref_file, action='r')
  do i = 1,size(var)


    allocate(Anew(lx2, lx3), Aref(lx2, lx3))
    call hnew%read(var(i), Anew)
    call href%read(var(i), Aref)

    if (.not.all(ieee_is_finite(Aref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // var(i)
    if (.not.all(ieee_is_finite(Anew))) error stop "NON-FINITE: " // file_name(new_file) // " " // var(i)

    if(all(isclose(Aref, Anew, real(rtol), real(atol)))) then
      if(P%debug) print '(A)', "OK: input:precip " // var(i) // " " // new_file
    else
      bad = bad + 1

      write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH:precip: " // file_name(new_file) // " " // var(i), &
      ' max diff:', maxval(abs(Aref - Anew)), &
      ' max & min ref:', maxval(Aref), minval(Aref), ' max & min new:', maxval(Anew), minval(Anew)
    endif

    deallocate(Anew, Aref)

  end do

  call hnew%close()
  call href%close()

  call dateinc(new%dtE0, ymd, UTsec)

  t = t + new%dtE0

end do

if (bad == 0) print '(a,1x,a)', "OK: Efield: ", new%E0dir

end function check_Efield


end submodule compare_in_h5
