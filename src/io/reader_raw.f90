submodule (reader) reader_raw

use, intrinsic :: iso_fortran_env, only : real64

implicit none (type, external)

contains


module procedure get_simsize2_raw
integer :: u

if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid size from file:', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) llon,llat
close(u)
end procedure get_simsize2_raw


module procedure get_simsize3_raw
!! note that these are sizes *including ghost cells*
integer :: u

if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid size from file:', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) lx1
read(u) lx2all
if (present(lx3all)) read(u) lx3all
close(u)
end procedure get_simsize3_raw


module procedure get_grid2_raw
integer :: u

if (debug) print '(A,/,A)', 'READ 2D (B-perp, B-perp) grid:', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) mlonp,mlatp
close(u)
end procedure get_grid2_raw


module procedure get_Efield_raw
integer :: u
real(real64) :: flagtmp

if (debug) print *, 'READ electric field data from:  ',path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
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
end procedure get_Efield_raw


module procedure get_precip_raw
integer :: u

if (debug) print *, 'READ precipitation data from: ', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) Qp,E0p
close(u)
end procedure get_precip_raw


module procedure get_neutral2_raw
integer :: u

if (debug) print *, 'READ neutral 2D data from: ', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) dnO,dnN2,dnO2,dvnrho,dvnz,dTn
close(u)
end procedure get_neutral2_raw


module procedure get_neutral3_raw
integer :: u

if (debug) print *, 'READ neutral 3D data from: ', path

open(newunit=u, file=path, status='old', form='unformatted', access='stream', action='read')
read(u) dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall
close(u)
end procedure get_neutral3_raw

end submodule reader_raw
