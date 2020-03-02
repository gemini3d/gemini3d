submodule (reader) reader_hdf5

use h5fortran, only: hdf5_file

implicit none

contains


module procedure get_simsize2
!! get x2 and x3 dimension sizes
type(hdf5_file) :: hf
character(:), allocatable :: fn
integer :: ierr
logical :: exists

if (index(path, 'simsize.h5') /= 0) then
  fn = path
else
  fn = path // '/simsize.h5'
endif
if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid size from file:', fn

inquire(file=fn, exist=exists)
if (.not.exists) then
   write(stderr,'(/,A,/,A,/)') 'ERROR: reader_hdf5:get_simsize2: generate grid with script--grid not present:',fn
   error stop 77
endif

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'get_simsize2: simsize.h5'

!> scripts can use variety of variable names
if(hf%exist('/llat', ierr)) then
  call hf%read('/llat', llat, ierr)
elseif(hf%exist('/Nlat', ierr)) then
  call hf%read('/Nlat', llat, ierr)
elseif(hf%exist('/lx2', ierr)) then
  call hf%read('/lx2', llat, ierr)
else
  error stop 'reader_hdf5:get_simsize2: llat / lx2'
endif
if(ierr/=0) error stop 'reader_hdf5:get_simsize2: llat / lx2'

if(hf%exist('/llon', ierr)) then
  call hf%read('/llon', llon, ierr)
elseif(hf%exist('/Nlat', ierr)) then
  call hf%read('/Nlat', llon, ierr)
elseif(hf%exist('/lx3', ierr)) then
  call hf%read('/lx3', llon, ierr)
else
  error stop 'reader_hdf5:get_simsize2: llon / lx3'
endif
if(ierr/=0) error stop 'reader_hdf5:get_simsize2: llon / lx3'

call hf%finalize(ierr)

end procedure get_simsize2


module procedure get_simsize3
!! get x1, x2, x3 dimension sizes
!! sizes include Ghost Cells
type(hdf5_file) :: hf
character(:), allocatable :: fn
integer :: ierr
logical :: exists

if (index(path, 'simsize.h5') /= 0) then
  fn = path
else
  fn = path // '/simsize.h5'
endif
if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid  size from file:', fn

inquire(file=fn, exist=exists)
if (.not.exists) then
   write(stderr,'(A,/,A)') 'ERROR: generate grid with script--grid not present: ',fn
   error stop 77
endif

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not open simsize.h5'

if (hf%exist("/lx1", ierr)) then
  call hf%read('/lx1', lx1, ierr)
  if(ierr/=0) error stop 'get_simsize3: lx1'
  call hf%read('/lx2', lx2all, ierr)
  if(ierr/=0) error stop 'get_simsize3: lx2'
  if (present(lx3all)) then
    call hf%read('/lx3', lx3all, ierr)
    if(ierr/=0) error stop 'get_simsize3: lx3'
  endif
elseif (hf%exist("/lxs", ierr)) then
  block
    integer :: lx(3)
    call hf%read("/lxs", lx, ierr)
    lx1 = lx(1)
    lx2all = lx(2)
    if (present(lx3all)) lx3all = lx(3)
  end block
elseif (hf%exist("/lx", ierr)) then
  block
    integer :: lx(3)
    call hf%read("/lx", lx, ierr)
    lx1 = lx(1)
    lx2all = lx(2)
    if (present(lx3all)) lx3all = lx(3)
  end block
endif

call hf%finalize(ierr)
if(ierr/=0) error stop 'could not close simsize.h5'

end procedure get_simsize3


module procedure get_grid2
type(hdf5_file) :: hf
character(:), allocatable :: fn
integer :: ierr

if (index(path, 'simgrid.h5') /= 0) then
  fn = path
else
  fn = path // '/simgrid.h5'
endif
if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid:', fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'get_grid2: could not open simgrid.h5'
call hf%read('/mlon', mlonp, ierr)
if(ierr/=0) error stop 'get_grid2: mlon'
call hf%read('/mlat', mlatp, ierr)
if(ierr/=0) error stop 'get_grid2: mlat'
call hf%finalize(ierr)

end procedure get_grid2


module procedure get_Efield
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn
integer :: ierr

fn = path // '.h5'
if (debug) print *, 'READ electric field data from file:  ',fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'get_Efield: could not read Efield HDF5'
call hf%read('/flagdirich', flagdirich, ierr)
if(ierr/=0) error stop 'get_Efield: flagdirich'
!! handle degenerate cases to avoid
!! "Operating system error: Cannot allocate memory Memory allocation failed"
if (size(E0xp, 1)==1) then
  call hf%read('/Exit', E0xp(1,:), ierr)
  if(ierr/=0) error stop
  call hf%read('/Eyit', E0yp(1,:), ierr)
  if(ierr/=0) error stop
  call hf%read('/Vminx1it', Vminx1p(1,:),ierr)
  if(ierr/=0) error stop
  call hf%read('/Vmaxx1it', Vmaxx1p(1,:), ierr)
  if(ierr/=0) error stop
elseif (size(E0xp, 2)==1) then
  call hf%read('/Exit', E0xp(:,1), ierr)
  if(ierr/=0) error stop
  call hf%read('/Eyit', E0yp(:,1), ierr)
  if(ierr/=0) error stop
  call hf%read('/Vminx1it', Vminx1p(:,1),ierr)
  if(ierr/=0) error stop
  call hf%read('/Vmaxx1it', Vmaxx1p(:,1), ierr)
  if(ierr/=0) error stop
else !< 3D
  call hf%read('/Exit', E0xp, ierr)
  if(ierr/=0) error stop
  call hf%read('/Eyit', E0yp, ierr)
  if(ierr/=0) error stop
  call hf%read('/Vminx1it', Vminx1p,ierr)
  if(ierr/=0) error stop
  call hf%read('/Vmaxx1it', Vmaxx1p, ierr)
  if(ierr/=0) error stop
endif
!! background fields and top/bottom boundary conditions
call hf%read('/Vminx2ist', Vminx2pslice,ierr)
if(ierr/=0) error stop
call hf%read('/Vmaxx2ist', Vmaxx2pslice, ierr)
if(ierr/=0) error stop
!! these only used for 3D simulations
call hf%read('/Vminx3ist', Vminx3pslice,ierr)
if(ierr/=0) error stop
call hf%read('/Vmaxx3ist', Vmaxx3pslice, ierr)
if(ierr/=0) error stop

call hf%finalize(ierr)
if(ierr/=0) error stop

end procedure get_Efield


module procedure get_precip
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn
integer :: ierr

fn = path // '.h5'
if (debug) print *, 'READ precipitation data from file:  ',fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not open precipitation HDF5'
!! handle degenerate cases to avoid
!! "Operating system error: Cannot allocate memory Memory allocation failed"
if (size(Qp, 1)==1) then
  call hf%read('/Qp', Qp(1,:), ierr)
  if(ierr/=0) error stop
  call hf%read('/E0p', E0p(1,:), ierr)
  if(ierr/=0) error stop
elseif (size(Qp, 2)==1) then
  call hf%read('/Qp', Qp(:,1), ierr)
  if(ierr/=0) error stop
  call hf%read('/E0p', E0p(:,1), ierr)
  if(ierr/=0) error stop
else  !< 3-D
  call hf%read('/Qp', Qp, ierr)
  if(ierr/=0) error stop
  call hf%read('/E0p', E0p, ierr)
  if(ierr/=0) error stop
endif

call hf%finalize(ierr)
if(ierr/=0) error stop

end procedure get_precip


module procedure get_neutral2
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn
integer :: ierr

fn = path // '.h5'
if (debug) print *, 'READ neutral 2D data from file:  ',fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'get_neutral2: could not open precipitation HDF5'

call hf%read('/dn0all', dnO, ierr)
if(ierr/=0) error stop 'get_neutral2: dn0'
call hf%read('/dnN2all', dnN2, ierr)
if(ierr/=0) error stop 'get_neutral2: dnN2'
call hf%read('/dnO2all', dnO2, ierr)
if(ierr/=0) error stop 'get_neutral2: dnO2'
call hf%read('/dvnrhoall', dvnrho, ierr)
if(ierr/=0) error stop 'get_neutral2: dvnrho'
call hf%read('/dvnzall', dvnz, ierr)
if(ierr/=0) error stop 'get_neutral2: dvnz'
call hf%read('/dTnall', dTn, ierr)
if(ierr/=0) error stop 'get_neutral2: dTn'

call hf%finalize(ierr)

end procedure get_neutral2


module procedure get_neutral3
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn
integer :: ierr

fn = path // '.h5'
if (debug) print *, 'READ neutral 3D data from file:  ',fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not open precipitation HDF5'

call hf%read('/dn0all', dnOall, ierr)
if(ierr/=0) error stop 'get_neutral3: dn0all'
call hf%read('/dnN2all', dnN2all, ierr)
if(ierr/=0) error stop 'get_neutral3: dnN2'
call hf%read('/dnO2all', dnO2all, ierr)
if(ierr/=0) error stop 'get_neutral3: dnO2all'
call hf%read('/dvnxall', dvnxall, ierr)
if(ierr/=0) error stop 'get_neutral3: dnvxall'
call hf%read('/dvnrhoall', dvnrhoall, ierr)
if(ierr/=0) error stop 'get_neutral3: dvnrhoall'
call hf%read('/dvnzall', dvnzall, ierr)
if(ierr/=0) error stop 'get_neutral3: dvnzall'
call hf%read('/dTnall', dTnall, ierr)
if(ierr/=0) error stop 'get_neutral3: dTnall'

call hf%finalize(ierr)

end procedure get_neutral3



end submodule reader_hdf5