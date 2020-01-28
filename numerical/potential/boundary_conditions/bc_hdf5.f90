submodule (potentialBCs_mumps) bc_raw

use h5fortran, only: hdf5_file

implicit none

contains


module procedure get_simsize
type(hdf5_file) :: hf
character(:), allocatable :: fn
integer :: ierr

fn = path // '/simsize.h5'
if (debug) print '(A,/,A)', 'Inputting electric field data size from file:', fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not read simsize.h5'
call hf%read('/llat', llat, ierr)
if(ierr/=0) error stop 'could not read llat'
call hf%read('/llon', llon, ierr)
if(ierr/=0) error stop 'could not read llon'
call hf%finalize(ierr)

end procedure get_simsize


module procedure get_simgrid
type(hdf5_file) :: hf
character(:), allocatable :: fn
integer :: ierr

fn = path // '/simgrid.hf'
if (debug) print '(A,/,A)', 'Inputting electric field grid:', fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not read simgrid.hf'
call hf%read('/mlon', mlonp, ierr)
if(ierr/=0) error stop 'could not read mlon'
call hf%read('/mlat', mlatp, ierr)
if(ierr/=0) error stop 'could not read mlat'
call hf%finalize(ierr)

end procedure get_simgrid


module procedure get_Efield
type(hdf5_file) :: hf
real(wp) :: flagtmp
character(:), allocatable :: fn
integer :: ierr

fn = path // '.h5'
if (debug) print *, 'Read: electric field data from file:  ',fn

call hf%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) error stop 'could not read Efield HDF5'
call hf%read('/flagdirich', flagdirich, ierr)
if(ierr/=0) error stop 'could not read flagdirich'
call hf%read('/Exit', E0xp, ierr)
if(ierr/=0) error stop 'could not read Ex'
call hf%read('/Eyit', E0yp, ierr)
if(ierr/=0) error stop 'could not read Ey'
call hf%read('/Vminx1it', Vminx1p,ierr)
call hf%read('/Vmaxx1it', Vmaxx1p, ierr)
!! background fields and top/bottom boundary conditions
call hf%read('/Vminx2ist', Vminx2pslice,ierr)
call hf%read('/Vmaxx2ist', Vmaxx2pslice, ierr)
!! these only used for 3D simulations
call hf%read('/Vminx3ist', Vminx3pslice,ierr)
call hf%read('/Vmaxx3ist', Vmaxx3pslice, ierr)
call hf%finalize(ierr)
end procedure get_Efield

end submodule bc_raw