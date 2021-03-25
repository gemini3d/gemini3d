module compare_h5
!! for safety, this trades efficiency for reliability
!! that is, we repeatedly open, close, allocate to help avoid
!! any weird bugs causing false positive/negative

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only : wp
use timeutils, only : date_filename, dateinc
use config, only : gemini_cfg, read_configfile
use h5fortran, only : hdf5_file
use pathlib, only : get_suffix, parent, file_name
use reader, only : get_simsize3
use assert, only : isclose

implicit none (type, external)

integer, parameter :: lsp=7

real(wp), parameter :: &
rtol = 1e-5_wp,  atol = 1e-8_wp, &
rtolJ = 0.01_wp, atolJ = 1e-7_wp, &
rtolV = 1e-5_wp, atolV = 50, &
rtolN = 1e-5_wp, atolN = 1e9_wp, &
rtolT = 1e-5_wp, atolT = 100

private
public :: check_plasma_hdf5

contains


subroutine check_plasma_hdf5(new_path, ref_path, all_ok)

type(gemini_cfg) :: cfg

character(*), intent(in) :: new_path, ref_path
logical, intent(out) :: all_ok

character(:), allocatable :: new_file, ref_file
integer :: i, ymd(3)

real(wp) :: UTsec
character(:), allocatable :: suffix
logical :: exists, ok

all_ok = .true.

cfg%infile = new_path // '/inputs/config.nml'
cfg%outdir = '.'
call read_configfile(cfg)

ymd = cfg%ymd0
UTsec = cfg%UTsec0
suffix = get_suffix(cfg%indatsize)

ref_file = date_filename(ref_path, ymd, UTsec) // suffix
inquire(file=ref_file, exist=exists)
if (.not. exists) then
  ! skip first file for old sim reference data
  call dateinc(cfg%dtout, ymd, UTsec)
  ref_file = date_filename(ref_path, ymd, UTsec) // suffix
  inquire(file=ref_file, exist=exists)
  if (.not. exists) error stop "reference data not found in " // ref_path
endif

do
  ref_file = date_filename(ref_path, ymd, UTsec) // suffix
  inquire(file=ref_file, exist=exists)
  if (.not. exists) exit
  !! last output file
  new_file = date_filename(new_path, ymd, UTsec) // suffix

  ok = checker(cfg, new_file, ref_file)
  if(.not. ok) then
    write(stderr,*) "gemini3d.compare: MISMATCHED data at", ymd, UTsec
    all_ok = .false.
  endif

  !! next time
  call dateinc(cfg%dtout, ymd, UTsec)
end do


end subroutine check_plasma_hdf5


logical function checker(cfg, new_file, ref_file)

type(gemini_cfg), intent(in) :: cfg
character(*), intent(in) :: new_file, ref_file

type(hdf5_file) :: hnew, href
integer :: i, bad
logical :: exists, ok

character(7), parameter :: varsT(2) = [character(7) :: 'Tavgall', 'TEall']
character(8), parameter :: varsV(3) = ['v1avgall', 'v2avgall', 'v3avgall']
character(5), parameter :: varsJ(3) = ['J1all', 'J2all', 'J3all']

real(wp) :: UThour1, UThour2, UTsec1, UTsec2
integer :: ymd1(3), ymd2(3)
integer :: lx1, lx2all, lx3all, Nlx1, Nlx2all, Nlx3all, flagoutput

call get_simsize3(parent(ref_file) // "/inputs", lx1, lx2all, lx3all)
call get_simsize3(parent(new_file) // "/inputs", Nlx1, Nlx2all, Nlx3all)

if(lx1 /= Nlx1) error stop 'lx1 not match ref: ' // new_file
if(lx2all /= Nlx2all) error stop 'lx2all not match ref: ' // new_file
if(lx3all /= Nlx3all) error stop 'lx3all not match ref: ' // new_file

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

flagoutput = -1
if (hnew%exist("/flagoutput")) then
  call hnew%read("/flagoutput", flagoutput)
elseif (href%exist("/flagoutput")) then
  call href%read("/flagoutput", flagoutput)
else
  flagoutput = cfg%flagoutput
endif

call href%read('/time/ymd', ymd1)
call hnew%read('/time/ymd', ymd2)

call href%read('/time/UThour', UThour1)
call hnew%read('/time/UThour', UThour2)

call hnew%finalize()
call href%finalize()

!> compare file simulation time
UTsec1 = UThour1*3600
UTsec2 = UThour2*3600
call dateinc(0._wp, ymd1, UTsec1)
call dateinc(0._wp, ymd2, UTsec2)
!! sanitize wrapping glitches due to non-integer timebase
!! due to non-integer timebase, can get hour-wrapping in Fortran code.
!! This would be fixed someday by using integer microsecond timebase

if (any(ymd1 /= ymd2)) error stop 'dates did not match: ' // new_file
if (abs(UTsec1 - UTsec2) > 0.1) error stop "UThour not match: " // new_file

bad = 0

if (cfg%flagoutput == 3) then
  !! just electron density
  bad = bad + check_var('neall', new_file, ref_file, rtolN, atolN, lx1, lx2all, lx3all)

elseif (cfg%flagoutput == 2) then

  bad = bad + check_var('neall', new_file, ref_file, rtolN, atolN, lx1, lx2all, lx3all)

  do i = 1,size(varsT)
    bad = bad + check_var(varsT(i), new_file, ref_file, rtolT, atolT, lx1, lx2all, lx3all)
  enddo

  do i = 1,size(varsV)
    bad = bad + check_var(varsV(i), new_file, ref_file, rtolV, atolV, lx1, lx2all, lx3all)
  enddo

  do i = 1,size(varsJ)
    bad = bad + check_var(varsJ(i), new_file, ref_file, rtolJ, atolJ, lx1, lx2all, lx3all)
  enddo

elseif(cfg%flagoutput==1) then

  do i = 1,size(varsJ)
    bad = bad + check_var(varsJ(i), new_file, ref_file, rtolJ, atolJ, lx1, lx2all, lx3all)
  enddo

  do i = 2,size(varsV)
    bad = bad + check_var(varsV(i), new_file, ref_file, rtolV, atolV, lx1, lx2all, lx3all)
  enddo

  !> Ne
  bad = bad + check_var('nsall', new_file, ref_file, rtolN, atolN, lx1, lx2all, lx3all, ionly=lsp, derived_name="ne")

  !> Te
  bad = bad + check_var('Tsall', new_file, ref_file, rtolT, atolT, lx1, lx2all, lx3all, ionly=lsp, derived_name="Te")

  !> Ti
  bad = bad + check_derived('Tsall', "Ti", new_file, ref_file, rtolT, atolT, lx1, lx2all, lx3all)

  !> v1
  bad = bad + check_derived('vs1all', "v1", new_file, ref_file, rtolV, atolV, lx1, lx2all, lx3all)

else
  error stop 'unknown flagoutput: ' // file_name(ref_file)
endif

checker = bad == 0

end function checker


integer function check_derived(name, derived_name, new_file, ref_file, rtol, atol, lx1, lx2all, lx3all) result(bad)

character(*), intent(in) :: new_file, ref_file
character(*), intent(in) :: name, derived_name
real(wp), intent(in) :: rtol, atol
integer, intent(in) :: lx1, lx2all, lx3all

real, dimension(:,:,:), allocatable :: D_new, D_ref
real, dimension(:,:,:,:), allocatable :: new, ref, ns_new, ns_ref

type(hdf5_file) :: hnew, href

bad = 0

allocate(new(lx1, lx2all, lx3all, lsp), ref(lx1, lx2all, lx3all, lsp))
allocate(ns_new(lx1, lx2all, lx3all, lsp), ns_ref(lx1, lx2all, lx3all, lsp))
allocate(D_ref(lx1, lx2all, lx3all), D_new(lx1, lx2all, lx3all))

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

call href%read('nsall', ns_ref)
call hnew%read('nsall', ns_new)
if (.not.all(ieee_is_finite(ns_ref))) error stop "NON-FINITE: " // file_name(ref_file) // " ns"
if (.not.all(ieee_is_finite(ns_new))) error stop "NON-FINITE: " // file_name(new_file) // " ns"

call hnew%read(name, new)
call href%read(name, ref)
if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(new_file) // " " // name
if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // name

call hnew%finalize()
call href%finalize()

D_ref = sum(ns_ref(:,:,:,1:6) * ref(:,:,:,1:6), dim=4) / ns_ref(:,:,:,LSP)
D_new = sum(ns_new(:,:,:,1:6) * new(:,:,:,1:6), dim=4) / ns_new(:,:,:,LSP)
if(all(isclose(D_ref, D_new, rtol, atol))) return

bad = bad + 1
write(stderr,*) "MISMATCH: " // file_name(new_file) // " ", derived_name, maxval(abs(D_ref - D_new))

end function check_derived


integer function check_var(name, new_file, ref_file, rtol, atol, lx1, lx2all, lx3all, ionly, derived_name) result(bad)

character(*), intent(in) :: new_file, ref_file
character(*), intent(in) :: name
real(wp), intent(in) :: rtol, atol
integer, intent(in) :: lx1, lx2all, lx3all
integer, intent(in), optional :: ionly
character(*), intent(in), optional :: derived_name

type(hdf5_file) :: hnew, href

real, dimension(:,:,:), allocatable :: new, ref
real, dimension(:,:,:,:), allocatable :: new4, ref4

if(present(ionly)) then
  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp))
endif

allocate(new(lx1, lx2all, lx3all), ref(lx1, lx2all, lx3all))


bad = 0

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

if(present(ionly)) then
  call hnew%read(name, new4)
  call href%read(name, ref4)

  new = new4(:,:,:,lsp)
  ref = ref4(:,:,:,lsp)
else
  call hnew%read(name, new)
  call href%read(name, ref)
endif

call hnew%finalize()
call href%finalize()

if (.not.all(ieee_is_finite(ref))) error stop "NON-FINITE: " // file_name(ref_file) // " " // name
if (.not.all(ieee_is_finite(new))) error stop "NON-FINITE: " // file_name(new_file) // " " // name

if(all(isclose(ref, new, rtol, atol))) return

bad = 1
if(present(derived_name)) then
  write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // derived_name, &
    ' max diff:', maxval(abs(ref - new)), ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(new)
else
  write(stderr,'(A,/,A,ES12.3,A,2ES12.3,A,2ES12.3)') "MISMATCH: " // file_name(new_file) // " " // name, &
    ' max diff:', maxval(abs(ref - new)), ' max & min ref:', maxval(ref), minval(ref), ' max & min new:', maxval(new), minval(ref)
endif



end function check_var


end module compare_h5
