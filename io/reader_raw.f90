submodule (reader) reader_raw

use, intrinsic :: iso_fortran_env, only : real64

implicit none

contains


module procedure get_simsize2
integer :: u
character(:), allocatable :: fn
logical :: exists

if (index(path, 'simsize.dat') /= 0) then
  fn = path
else
  fn = path // '/simsize.dat'
endif
if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid size from file:', fn

inquire(file=fn, exist=exists)
if (.not.exists) then
   write(stderr,'(A,/,A)') 'ERROR: generate grid with script--grid not present: ',fn
   error stop 77
endif

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) llon,llat
close(u)
end procedure get_simsize2


module procedure get_simsize3
!! note that these are sizes *including ghost cells*
integer :: u
character(:), allocatable :: fn
logical :: exists

if (index(path, 'simsize.dat') /= 0) then
  fn = path
else
  fn = path // '/simsize.dat'
endif
if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid size from file:', fn

inquire(file=fn, exist=exists)
if (.not.exists) then
   write(stderr,'(A,/,A)') 'ERROR: generate grid with script--grid not present: ',fn
   error stop 77
endif

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) lx1
read(u) lx2all
if (present(lx3all)) read(u) lx3all
close(u)
end procedure get_simsize3


module procedure get_grid2
integer :: u
character(:), allocatable :: fn

if (index(path, 'simgrid.dat') /= 0) then
  fn = path
else
  fn = path // '/simgrid.dat'
endif
if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid:', fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) mlonp,mlatp
close(u)
end procedure get_grid2


module procedure get_Efield
integer :: u
real(real64) :: flagtmp
character(:), allocatable :: fn

fn = path // '.dat'
if (debug) print *, 'READ electric field data from:  ',fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) flagtmp
!! NOTE: this is mistakenly a float from Matlab
!! to keep compatibility with old files, we left it as real64.
!! New work should be using HDF5 instead of raw in any case.
flagdirich = int(flagtmp,4)
read(u) E0xp,E0yp
read(u) Vminx1p,Vmaxx1p
!! background fields and top/bottom boundary conditions
read(u) Vminx2pslice,Vmaxx2pslice
!! these only used for 3D simulations
read(u) Vminx3pslice,Vmaxx3pslice
close(u)
end procedure get_Efield


module procedure get_precip
integer :: u
character(:), allocatable :: fn

fn = path // '.dat'
if (debug) print *, 'READ precipitation data from:  ',fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) Qp,E0p
close(u)
end procedure get_precip


module procedure get_neutral2
integer :: u
character(:), allocatable :: fn

fn = path // '.dat'
if (debug) print *, 'READ neutral 2D data from:  ',fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) dnO,dnN2,dnO2,dvnrho,dvnz,dTn
close(u)
end procedure get_neutral2


module procedure get_neutral3
integer :: u
character(:), allocatable :: fn

fn = path // '.dat'
if (debug) print *, 'READ neutral 3D data from:  ',fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
close(u)
end procedure get_neutral3

end submodule reader_raw