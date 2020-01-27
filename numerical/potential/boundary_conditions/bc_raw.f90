submodule (potentialBCs_mumps) bc_raw

implicit none

contains


module procedure get_simsize
integer :: u
character(:), allocatable :: fn

fn = path // '/simsize.dat'
if (debug) print '(A,/,A)', 'Inputting electric field data size from file:', fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) llon,llat
close(u)
end procedure get_simsize


module procedure get_simgrid
integer :: u
character(:), allocatable :: fn

fn = path // '/simgrid.dat'
if (debug) print '(A,/,A)', 'Inputting electric field grid:', fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) mlonp,mlatp
close(u)
end procedure get_simgrid


module procedure get_Efield
integer :: u
real(wp) :: flagtmp
character(:), allocatable :: fn

fn = path // '.dat'
if (debug) print *, 'Read: electric field data from file:  ',fn

open(newunit=u, file=fn, status='old', form='unformatted', access='stream')
read(u) flagtmp  !< FIXME: this is mistakenly a float from Matlab
flagdirich = int(flagtmp,4)
read(u) E0xp,E0yp
read(u) Vminx1p,Vmaxx1p
!! background fields and top/bottom boundary conditions
read(u) Vminx2pslice,Vmaxx2pslice
!! these only used for 3D simulations
read(u) Vminx3pslice,Vmaxx3pslice
close(u)
end procedure get_Efield

end submodule bc_raw