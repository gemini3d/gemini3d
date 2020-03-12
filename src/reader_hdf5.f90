submodule (reader) reader_hdf5

use h5fortran, only: hdf5_file

implicit none

contains


module procedure get_simsize2
!! get x2 and x3 dimension sizes
type(hdf5_file) :: hf
character(:), allocatable :: fn
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

call hf%initialize(fn, status='old',action='r')

!> scripts can use variety of variable names
if(hf%exist('/llat')) then
  call hf%read('/llat', llat)
elseif(hf%exist('/Nlat')) then
  call hf%read('/Nlat', llat)
elseif(hf%exist('/lx2')) then
  call hf%read('/lx2', llat)
else
  error stop 'reader_hdf5:get_simsize2: llat / lx2'
endif

if(hf%exist('/llon')) then
  call hf%read('/llon', llon)
elseif(hf%exist('/Nlat')) then
  call hf%read('/Nlat', llon)
elseif(hf%exist('/lx3')) then
  call hf%read('/lx3', llon)
else
  error stop 'reader_hdf5:get_simsize2: llon / lx3'
endif

call hf%finalize()

end procedure get_simsize2


module procedure get_simsize3
!! get x1, x2, x3 dimension sizes
!! sizes include Ghost Cells
type(hdf5_file) :: hf
character(:), allocatable :: fn
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

call hf%initialize(fn, status='old',action='r')

if (hf%exist("/lx1")) then
  call hf%read('/lx1', lx1)
  call hf%read('/lx2', lx2all)
  if (present(lx3all)) then
    call hf%read('/lx3', lx3all)
  endif
elseif (hf%exist("/lxs")) then
  block
    integer :: lx(3)
    call hf%read("/lxs", lx)
    lx1 = lx(1)
    lx2all = lx(2)
    if (present(lx3all)) lx3all = lx(3)
  end block
elseif (hf%exist("/lx")) then
  block
    integer :: lx(3)
    call hf%read("/lx", lx)
    lx1 = lx(1)
    lx2all = lx(2)
    if (present(lx3all)) lx3all = lx(3)
  end block
endif

call hf%finalize()

end procedure get_simsize3


module procedure get_grid2
type(hdf5_file) :: hf
character(:), allocatable :: fn

if (index(path, 'simgrid.h5') /= 0) then
  fn = path
else
  fn = path // '/simgrid.h5'
endif
if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid:', fn

call hf%initialize(fn, status='old',action='r')
call hf%read('/mlon', mlonp)
call hf%read('/mlat', mlatp)
call hf%finalize()

end procedure get_grid2


module procedure get_Efield
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn

fn = path // '.h5'
if (debug) print *, 'READ electric field data from file:  ',fn

call hf%initialize(fn, status='old',action='r')
call hf%read('/flagdirich', flagdirich)
!! TODO: fix Matlab to include degenerate dimension, to avoid
!! "Operating system error: Cannot allocate memory Memory allocation failed"
! if (size(E0xp, 1)==1) then
!   call hf%read('/Exit',  E0xp) ! E0xp(1,:))
!   call hf%read('/Eyit', E0yp(1,:))
!   call hf%read('/Vminx1it', Vminx1p(1,:))
!   call hf%read('/Vmaxx1it', Vmaxx1p(1,:))
! elseif (size(E0xp, 2)==1) then
!   call hf%read('/Exit', E0xp(:,1))
!   call hf%read('/Eyit', E0yp(:,1))
!   call hf%read('/Vminx1it', Vminx1p(:,1))
!   call hf%read('/Vmaxx1it', Vmaxx1p(:,1))
! else !< 3D
  call hf%read('/Exit', E0xp)
  call hf%read('/Eyit', E0yp)
  call hf%read('/Vminx1it', Vminx1p)
  call hf%read('/Vmaxx1it', Vmaxx1p)
! endif
!! background fields and top/bottom boundary conditions
call hf%read('/Vminx2ist', Vminx2pslice)
call hf%read('/Vmaxx2ist', Vmaxx2pslice)
!! these only used for 3D simulations
call hf%read('/Vminx3ist', Vminx3pslice)
call hf%read('/Vmaxx3ist', Vmaxx3pslice)

call hf%finalize()

end procedure get_Efield


module procedure get_precip
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn

fn = path // '.h5'
if (debug) print *, 'READ precipitation data from file:  ',fn

call hf%initialize(fn, status='old',action='r')
!! TODO: fix Matlab to include degenerate dimension, to avoid
!! "Operating system error: Cannot allocate memory Memory allocation failed"
! if (size(Qp, 1)==1) then
!   call hf%read('/Qp', Qp(1,:))
!   call hf%read('/E0p', E0p(1,:))
! elseif (size(Qp, 2)==1) then
!   call hf%read('/Qp', Qp(:,1))
!   call hf%read('/E0p', E0p(:,1))
! else  !< 3-D
  call hf%read('/Qp', Qp)
  call hf%read('/E0p', E0p)
! endif

call hf%finalize()

end procedure get_precip


module procedure get_neutral2
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn

fn = path // '.h5'
if (debug) print *, 'READ neutral 2D data from file:  ',fn

call hf%initialize(fn, status='old',action='r')

call hf%read('/dn0all', dnO)
call hf%read('/dnN2all', dnN2)
call hf%read('/dnO2all', dnO2)
call hf%read('/dvnrhoall', dvnrho)
call hf%read('/dvnzall', dvnz)
call hf%read('/dTnall', dTn)

call hf%finalize()

end procedure get_neutral2


module procedure get_neutral3
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn

fn = path // '.h5'
if (debug) print *, 'READ neutral 3D data from file:  ',fn

call hf%initialize(fn, status='old',action='r')

call hf%read('/dn0all', dnOall)
call hf%read('/dnN2all', dnN2all)
call hf%read('/dnO2all', dnO2all)
call hf%read('/dvnxall', dvnxall)
call hf%read('/dvnrhoall', dvnrhoall)
call hf%read('/dvnzall', dvnzall)
call hf%read('/dTnall', dTnall)

call hf%finalize()

end procedure get_neutral3



end submodule reader_hdf5