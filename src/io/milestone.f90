submodule (io) milestone

use timeutils, only : date_filename,dateinc
use h5fortran, only : h5exist, hdf5_file
use hdf5, only: H5T_NATIVE_DOUBLE, h5tequal_f
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

implicit none (type,external)   !! external procedures must be explicitly denoted thusly

contains

module procedure find_milestone

!> search path having output rate cadence (s) and find the last file that is a milestone.
integer, dimension(3) :: ymd
real(wp) :: UTsec
character(:), allocatable :: fn
logical :: exists
logical :: is_double
real(wp) :: tsim, saved_ut
integer :: saved_ymd(3),schema,complete,realbits,i
integer :: type_error
character(3), parameter :: core_fields(4)=[character(3)::'ns','Ts','vs1','Phi']
type(hdf5_file) :: hf

tsim = 0
tmile = 0

ymd = cfg%ymd0
UTsec = cfg%UTsec0
ymdmile = cfg%ymd0
UTsecmile = cfg%UTsec0

filemile = date_filename(cfg%outdir, ymd, UTsec) // ".h5"
!! This presumes the first file output is a milestone.
!! We don't test the situation wheere a first output was not produced.
!! User should not be restarting in that case.

if (cfg%mcadence <= 0 .and. cfg%flagoutput/=1) then      !okay for milestone if full output specified
!! milestone was not in config.nml
  inquire(file=filemile, exist=exists)
 ! error stop filemile
  if (exists) error stop 'a fresh simulation should not have data in output directory: ' // filemile
  return
endif

milesearch : do
  !! new filename, add the 1 if it is the first
  fn = date_filename(cfg%outdir, ymd, UTsec) // ".h5"

  inquire(file=fn, exist=exists)
  if ( .not. exists ) exit milesearch
  !! last output file

  if (h5exist(fn, '/nsall')) then
    ! A partially written frame is not a restart checkpoint.
    if (.not.h5exist(fn,'/Tsall') .or. .not.h5exist(fn,'/vs1all') .or. &
        .not.h5exist(fn,'/Phiall') .or. .not.h5exist(fn,'/time/ymd') .or. &
        .not.h5exist(fn,'/time/UThour')) error stop 'Incomplete restart frame: ' // fn
    if (h5exist(fn,'/restart_core')) then
      if (.not.h5exist(fn,'/restart_core/complete') .or. .not.h5exist(fn,'/restart_core/schema') .or. &
          .not.h5exist(fn,'/restart_core/realbits') .or. .not.h5exist(fn,'/restart_core/ns') .or. &
          .not.h5exist(fn,'/restart_core/Ts') .or. .not.h5exist(fn,'/restart_core/vs1') .or. &
          .not.h5exist(fn,'/restart_core/Phi')) error stop 'Incomplete core restart record: '//fn
      call hf%open(fn,action='r')
      call hf%read('/restart_core/schema',schema)
      call hf%read('/restart_core/complete',complete)
      call hf%read('/restart_core/realbits',realbits)
      if (schema/=1.or.complete/=1.or.realbits/=storage_size(1._wp)) &
        error stop 'Unsupported core restart schema or precision: '//fn
      if (hf%ndim('/restart_core/Phi')/=3) error stop 'Core restart potential must be full 3D: '//fn
      do i=1,size(core_fields)
        call h5tequal_f(hf%dtype('/restart_core/'//trim(core_fields(i))),H5T_NATIVE_DOUBLE,is_double,type_error)
        if (type_error/=0) error stop 'Cannot compare core restart datatype: '//fn
        if (.not.is_double) error stop 'Core restart fields must contain float64 data: '//fn
      enddo
      call hf%close()
    elseif (cfg%potsolve==3) then
      error stop 'Field-resolved restart requires a full 3D potential checkpoint; saved slab is insufficient.'
    endif
    call hf%open(fn,action='r')
    call hf%read('/time/ymd',saved_ymd)
    call hf%read('/time/UThour',saved_ut)
    call hf%close()
    saved_ut=saved_ut*3600._wp
    if (.not.ieee_is_finite(saved_ut)) error stop 'Nonfinite restart time: ' // fn
    if (any(saved_ymd/=ymd) .or. abs(saved_ut-UTsec)>0.005_wp) error stop 'Restart timestamp mismatch: ' // fn
    !! this file is milestone
    ymdmile=ymd
    UTsecmile=UTsec
    filemile=fn
    tmile=tsim
  end if

  !! next time
  call dateinc(cfg%dtout, ymd,UTsec)
  tsim = tsim + cfg%dtout
end do milesearch

end procedure find_milestone

end submodule milestone
