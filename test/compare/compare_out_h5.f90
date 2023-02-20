submodule (compare_h5) compare_out_h5

implicit none (type, external)

contains

module procedure check_plasma_output_hdf5

type(gemini_cfg) :: cfg

integer :: ymd(3), lx1, lx2all, lx3all

real(wp) :: UTsec, t
logical :: ok

call check_simsize(new_path, ref_path, lx1, lx2all, lx3all)

check_plasma_output_hdf5 = .true.

cfg%infile = new_path // '/inputs/config.nml'
cfg%outdir = '.'  !< not used, just to pass checks
call read_configfile(cfg)

ymd = cfg%ymd0
UTsec = cfg%UTsec0
t = 0

do while (t <= cfg%tdur)

  ok = check_out(cfg, &
    new_file=date_filename(new_path, ymd, UTsec) // suffix(cfg%indatsize), &
    ref_file=date_filename(ref_path, ymd, UTsec) // suffix(cfg%indatsize), &
    lx1=lx1, lx2all=lx2all, lx3all=lx3all, P=P)

  if(.not. ok) then
    write(stderr,*) "gemini3d.compare: MISMATCHED data at", ymd, UTsec
    check_plasma_output_hdf5 = .false.
  endif

  !! next time
  call dateinc(cfg%dtout, ymd, UTsec)
  t = t + cfg%dtout
end do

end procedure check_plasma_output_hdf5


logical function check_out(cfg, new_file, ref_file, lx1, lx2all, lx3all, P)

type(gemini_cfg), intent(in) :: cfg
character(*), intent(in) :: new_file, ref_file
integer, intent(in) :: lx1, lx2all, lx3all
class(params), intent(in) :: P

type(hdf5_file) :: hnew, href
integer :: i, bad

character(7), parameter :: varsT(2) = [character(7) :: 'Tavgall', 'TEall']
character(8), parameter :: varsV(3) = ['v1avgall', 'v2avgall', 'v3avgall']
character(5), parameter :: varsJ(3) = ['J1all', 'J2all', 'J3all']

integer :: flagoutput

if(.not. is_file(new_file)) error stop "gemini3d.compare:check_out: new data file not found: " // new_file
if(.not. is_file(ref_file)) error stop "gemini3d.compare:check_out: reference data file not found: " // ref_file

call hnew%open(new_file, action='r')
call href%open(ref_file, action='r')

if (hnew%exist("/flagoutput")) then
  call hnew%read("/flagoutput", flagoutput)
elseif (href%exist("/flagoutput")) then
  call href%read("/flagoutput", flagoutput)
else
  flagoutput = cfg%flagoutput
endif

call check_time(hnew, href)

bad = 0

select case (flagoutput)
case (3)
  !! just electron density
  bad = bad + check_var('neall', hnew, href, rtolN, atolN, lx1, lx2all, lx3all, P)

case (2)

  bad = bad + check_var('neall', hnew, href, rtolN, atolN, lx1, lx2all, lx3all, P)

  do i = 1,size(varsT)
    bad = bad + check_var(varsT(i), hnew, href, rtolT, atolT, lx1, lx2all, lx3all, P)
  enddo

  do i = 1,size(varsV)
    bad = bad + check_var(varsV(i), hnew, href, rtolV, atolV, lx1, lx2all, lx3all, P)
  enddo

  do i = 1,size(varsJ)
    bad = bad + check_var(varsJ(i), hnew, href, rtolJ, atolJ, lx1, lx2all, lx3all, P)
  enddo

case (1)

  do i = 1,size(varsJ)
    bad = bad + check_var(varsJ(i), hnew, href, rtolJ, atolJ, lx1, lx2all, lx3all, P)
  enddo

  do i = 2,size(varsV)
    bad = bad + check_var(varsV(i), hnew, href, rtolV, atolV, lx1, lx2all, lx3all, P)
  enddo

  !> Ne
  bad = bad + check_var('nsall', hnew, href, rtolN, atolN, lx1, lx2all, lx3all, P, ionly=lsp, derived_name="ne")

  !> Te
  bad = bad + check_var('Tsall', hnew, href, rtolT, atolT, lx1, lx2all, lx3all, P, ionly=lsp, derived_name="Te")

  !> Ti
  bad = bad + check_derived('Tsall', "Ti", hnew, href, rtolT, atolT, lx1, lx2all, lx3all, P)

  !> v1
  bad = bad + check_derived('vs1all', "v1", hnew, href, rtolV, atolV, lx1, lx2all, lx3all, P)

case default
  error stop 'unknown flagoutput: ' // file_name(href%filename)
end select

check_out = bad == 0

end function check_out


integer function check_derived(name, derived_name, hnew, href, rtol, atol, lx1, lx2all, lx3all, P) result(bad)

type(hdf5_file), intent(in) :: hnew, href
character(*), intent(in) :: name, derived_name
real(wp), intent(in) :: rtol, atol
integer, intent(in) :: lx1, lx2all, lx3all
class(params), intent(in) :: P

real, dimension(:,:,:), allocatable :: D_new, D_ref
real, dimension(:,:,:,:), allocatable :: new, ref, ns_new, ns_ref

bad = 0

allocate(new(lx1, lx2all, lx3all, lsp), ref(lx1, lx2all, lx3all, lsp))
allocate(ns_new(lx1, lx2all, lx3all, lsp), ns_ref(lx1, lx2all, lx3all, lsp))
allocate(D_ref(lx1, lx2all, lx3all), D_new(lx1, lx2all, lx3all))

call href%read('nsall', ns_ref)
call hnew%read('nsall', ns_new)
if (.not.all(ieee_is_finite(ns_ref))) error stop "NON-FINITE: " // file_name(href%filename) // " ns"
if (.not.all(ieee_is_finite(ns_new))) error stop "NON-FINITE: " // file_name(hnew%filename) // " ns"

call hnew%read(name, new)
call href%read(name, ref)
if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(hnew%filename) // " " // name
if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(href%filename) // " " // name

D_ref = sum(ns_ref(:,:,:,1:6) * ref(:,:,:,1:6), dim=4) / ns_ref(:,:,:,LSP)
D_new = sum(ns_new(:,:,:,1:6) * new(:,:,:,1:6), dim=4) / ns_new(:,:,:,LSP)

if(all(isclose(D_ref, D_new, real(rtol), real(atol)))) then
  if(P%debug) print '(A)', "OK: output: " // derived_name // " " // hnew%filename
  return
endif

bad = 1
write(stderr,*) "MISMATCH: " // file_name(hnew%filename) // " ", derived_name, maxval(abs(D_ref - D_new))

call plot_diff(hnew%filename, href%filename, derived_name, "out", P)

end function check_derived


integer function check_var(name, hnew, href, rtol, atol, lx1, lx2all, lx3all, P, ionly, derived_name) result(bad)

type(hdf5_file), intent(in) :: hnew, href
character(*), intent(in) :: name
real(wp), intent(in) :: rtol, atol
integer, intent(in) :: lx1, lx2all, lx3all
class(params), intent(in) :: P
integer, intent(in), optional :: ionly
character(*), intent(in), optional :: derived_name

real, dimension(:,:,:), allocatable :: new, ref
real, dimension(:,:,:,:), allocatable :: new4, ref4

if(present(ionly)) then
  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp))
endif

allocate(new(lx1, lx2all, lx3all), ref(lx1, lx2all, lx3all))


bad = 0

if(present(ionly)) then
  call hnew%read(name, new4)
  call href%read(name, ref4)

  new = new4(:,:,:,lsp)
  ref = ref4(:,:,:,lsp)
else
  call hnew%read(name, new)
  call href%read(name, ref)
endif

if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(href%filename) // " " // name
if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(hnew%filename) // " " // name

if(all(isclose(ref, new, real(rtol), real(atol)))) then
  if(P%debug) then
    if (present(derived_name)) then
      print '(A)', "OK: output: " // derived_name
    else
      print '(A)', "OK: output: " // name
    endif
  endif
  return
endif

!> mismatch message
bad = 1
if(present(derived_name)) then
  write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(hnew%filename) // " " // derived_name, &
    ' max diff:', maxval(abs(ref - new)), ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(new)
else
  write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(hnew%filename) // " " // name, &
    ' max diff:', maxval(abs(ref - new)), ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(new)
endif

!> optional plotting
call plot_diff(hnew%filename, href%filename, name, "out", P)


end function check_var

end submodule compare_out_h5
