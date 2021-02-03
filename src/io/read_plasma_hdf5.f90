module read_plasma_h5

use phys_consts, only : wp
use timeutils, only : date_filename, dateinc
use config, only : gemini_cfg, read_configfile
use h5fortran, only : hdf5_file
use pathlib, only : get_suffix, parent, file_name
use reader, only : get_simsize3
use assert, only : assert_allclose

implicit none (type, external)

contains


subroutine read_plasma_hdf5(new_path, ref_path)

type(hdf5_file) :: hnew, href
type(gemini_cfg) :: cfg

character(*), intent(in) :: new_path, ref_path
character(:), allocatable :: new_file, ref_file
integer :: i, ymd(3)

real(wp) :: UTsec
character(:), allocatable :: suffix
logical :: exists

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

  call checker(cfg, new_file, ref_file)

  !! next time
  call dateinc(cfg%dtout, ymd, UTsec)
end do


end subroutine read_plasma_hdf5


subroutine checker(cfg, new_file, ref_file)

type(gemini_cfg), intent(in) :: cfg
character(*), intent(in) :: new_file, ref_file

type(hdf5_file) :: hnew, href
integer :: i
logical :: exists

real, allocatable :: new2(:,:), ref2(:,:)
real, dimension(:,:,:), allocatable :: new3, ref3, Ti_new, Ti_ref, v1_new, v1_ref
real, dimension(:,:,:,:), allocatable :: new4, ref4, ns_new, ns_ref

character(7), parameter :: varsT(2) = [character(7) :: 'Tavgall', 'TEall']
character(8), parameter :: varsV(3) = ['v1avgall', 'v2avgall', 'v3avgall']
character(5), parameter :: varsJ(3) = ['J1all', 'J2all', 'J3all']

character(:), allocatable :: k, rn

real(wp) :: UThour1, UThour2
integer :: ymd1(3), ymd2(3)
integer :: lx1, lx2all, lx3all, Nlx1, Nlx2all, Nlx3all, lsp

real(wp), parameter :: &
rtol = 1e-5_wp,     atol = 1e-8_wp, &
rtolJ = 1e-5_wp, atolJ = 1e-7_wp, &
rtolV = 1e-5_wp, atolV = 50, &
rtolN = 1e-5_wp, atolN = 1e9_wp, &
rtolT = 1e-5_wp, atolT = 100

lsp = 7

rn = file_name(ref_file)

call get_simsize3(parent(ref_file) // "/inputs", lx1, lx2all, lx3all)
call get_simsize3(parent(new_file) // "/inputs", Nlx1, Nlx2all, Nlx3all)

if(lx1 /= Nlx1) error stop 'lx1 not match ref: ' // new_file
if(lx2all /= Nlx2all) error stop 'lx2all not match ref: ' // new_file
if(lx3all /= Nlx3all) error stop 'lx3all not match ref: ' // new_file

if(lx3all /= 1) then
  !! not swapped
  allocate(new4(lx1, lx2all, lx3all, lsp), ref4(lx1, lx2all, lx3all, lsp), &
    ns_new(lx1, lx2all, lx3all, lsp), ns_ref(lx1, lx2all, lx3all, lsp))
  allocate(new3(lx1, lx2all, lx3all), ref3(lx1, lx2all, lx3all), &
    Ti_ref(lx1, lx2all, lx3all), Ti_new(lx1, lx2all, lx3all), v1_ref(lx1, lx2all, lx3all), v1_new(lx1, lx2all, lx3all))
else
  !! 2D sim, swapped
  allocate(new4(lx1, lx3all, lx2all, lsp), ref4(lx1, lx3all, lx2all, lsp), &
    ns_new(lx1, lx3all, lx2all, lsp), ns_ref(lx1, lx3all, lx2all, lsp))
  allocate(new3(lx1, lx3all, lx2all), ref3(lx1, lx3all, lx2all), &
     Ti_ref(lx1, lx3all, lx2all), Ti_new(lx1, lx3all, lx2all), v1_ref(lx1, lx3all, lx2all), v1_new(lx1, lx3all, lx2all))
endif

call hnew%initialize(new_file, status='old',action='r')
call href%initialize(ref_file, status='old',action='r')

call href%read('/time/ymd', ymd1)
call hnew%read('/time/ymd', ymd2)

if (any(ymd1 /= ymd2)) error stop 'dates did not match: ' // new_file

call href%read('/time/UThour', UThour1)
call hnew%read('/time/UThour', UThour2)

if (abs(UThour1 - UThour2) > 0.01) error stop "UThour not match: " // new_file

! print *, "flagoutput:", cfg%flagoutput

if (cfg%flagoutput == 3) then
  !! just electron density
  k = 'neall'
  call hnew%read(k, new3)
  call href%read(k, ref3)
  call assert_allclose(ref3, new3, rtolN, atolN, err_msg= rn // " " // k)

elseif (cfg%flagoutput == 2) then

  k = 'neall'
  call hnew%read(k, new3)
  call href%read(k, ref3)
  call assert_allclose(ref3, new3, rtolN, atolN, err_msg= rn // " " // k)

  do i = 1,size(varsT)
    k = varsT(i)
    call hnew%read(k, new3)
    call href%read(k, ref3)
    call assert_allclose(ref3, new3, rtolT, atolT, err_msg= rn // " " // k)
  enddo

  do i = 1,size(varsV)
    k = varsV(i)
    call hnew%read(k, new3)
    call href%read(k, ref3)
    call assert_allclose(ref3, new3, rtolV, atolV, err_msg= rn // " " // k)
  enddo

  do i = 1,size(varsJ)
    k = varsJ(i)
    call hnew%read(k, new3)
    call href%read(k, ref3)
    call assert_allclose(ref3, new3, rtolJ, atolJ, err_msg= rn // " " // k)
  enddo

elseif(cfg%flagoutput==1) then

  do i = 1,size(varsJ)
    k = varsJ(i)
    call hnew%read(k, new3)
    call href%read(k, ref3)
    call assert_allclose(ref3, new3, rtolJ, atolJ, err_msg= rn // " " // k)
  enddo

  do i = 2,size(varsV)
    k = varsV(i)
    call hnew%read(k, new3)
    call href%read(k, ref3)
    call assert_allclose(ref3, new3, rtolV, atolV, err_msg= rn // " " // k)
  enddo
  k = 'vs1all'
  call hnew%read(k, new4)
  call href%read(k, ref4)
  !> v1
  v1_ref = sum(ns_ref(:,:,:,1:6) * ref4(:,:,:,1:6), dim=4) / ns_ref(:,:,:,LSP)
  v1_new = sum(ns_new(:,:,:,1:6) * new4(:,:,:,1:6), dim=4) / ns_new(:,:,:,LSP)
  call assert_allclose(ref4, new4, rtolV, atolV, err_msg= rn // " " // k)

  k = 'nsall'
  call hnew%read(k, ns_new)
  call href%read(k, ns_ref)
  !> Ne
  call assert_allclose(ns_ref(:,:,:,LSP), ns_new(:,:,:,LSP), rtolN, atolN, err_msg= rn // " Ne")

  k = 'Tsall'
  call hnew%read(k, new4)
  call href%read(k, ref4)
  !> Te
  call assert_allclose(ref4(:,:,:,lsp), new4(:,:,:,lsp), rtolT, atolT, err_msg= rn // " Te")
  !> Ti
  Ti_ref = sum(ns_ref(:,:,:,1:6) * ref4(:,:,:,1:6), dim=4) / ns_ref(:,:,:,LSP)
  Ti_new = sum(ns_new(:,:,:,1:6) * new4(:,:,:,1:6), dim=4) / ns_new(:,:,:,LSP)
  call assert_allclose(Ti_ref, Ti_new, rtolT, atolT, err_msg= rn // " Ti")

else
  error stop 'unknown flagoutput: ' // rn
endif

call hnew%finalize()
call href%finalize()

end subroutine checker


end module read_plasma_h5
