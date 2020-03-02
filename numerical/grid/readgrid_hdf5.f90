submodule (grid:grid_read) readgrid_hdf5

use phys_consts, only: debug
use h5fortran, only: hdf5_file

implicit none

contains

module procedure get_grid3

type(hdf5_file) :: hout
character(:), allocatable :: fn
integer :: ierr

if (index(path, 'simgrid.h5') /= 0) then
  fn = path
else
  fn = path // '/simgrid.h5'
endif
if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid:', fn

call hout%initialize(fn, ierr, status='old',action='r')
if(ierr/=0) then
  write(stderr,'(A,/,A)') 'get_grid3: could not open simgrid HDF5 file ', fn
  error stop
endif

!> reads common to 2D and 3D
call hout%read('/x1', x%x1, ierr)
if(ierr/=0) error stop
call hout%read('/x1i', x%x1i, ierr)
if(ierr/=0) error stop
call hout%read('/dx1b', x%dx1, ierr)
if(ierr/=0) error stop
call hout%read('/dx1h', x%dx1i, ierr)
if(ierr/=0) error stop

if (flagswap/=1) then
  !! normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension

  call hout%read('/x2', x%x2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/x2i', x%x2iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx2b', x%dx2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx2h', x%dx2iall, ierr)
  if(ierr/=0) error stop

  call hout%read('/x3', x%x3all, ierr)
  if(ierr/=0) error stop
  call hout%read('/x3i', x%x3iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx3b', x%dx3all, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx3h', x%dx3iall, ierr)
  if(ierr/=0) error stop

  call hout%read('/h1', x%h1all, ierr)
  if(ierr/=0) error stop
  call hout%read('/h2', x%h2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/h3', x%h3all, ierr)
  if(ierr/=0) error stop

  call hout%read('/h1x1i', x%h1x1iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h2x1i', x%h2x1iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h3x1i', x%h3x1iall, ierr)
  if(ierr/=0) error stop

  call hout%read('/h1x2i', x%h1x2iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h2x2i', x%h2x2iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h3x2i', x%h3x2iall, ierr)
  if(ierr/=0) error stop

  call hout%read('/h1x3i', x%h1x3iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h2x3i', x%h2x3iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/h3x3i', x%h3x3iall, ierr)
  if(ierr/=0) error stop

  call hout%read('/gx1', g1all, ierr)
  if(ierr/=0) error stop
  call hout%read('/gx2', g2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/gx3', g3all, ierr)
  if(ierr/=0) error stop

  call hout%read('/alt', altall, ierr)
  if(ierr/=0) error stop
  call hout%read('/glat', glatall, ierr)
  if(ierr/=0) error stop
  call hout%read('/glon', glonall, ierr)
  if(ierr/=0) error stop

  call hout%read('/Bmag', Bmagall, ierr)
  if(ierr/=0) error stop
  !! corner case: for 2d sims, Incall has degenerate dimension
  !! and HDF5 file may be a 1-D vector, which can cause glibc
  !! error "corrupted size vs. prev_size" and hard crash.
  !! Since it's just one variable, easier to check and handle
  if (lx2all==1) then
    call hout%read('/I', Incall(1,:), ierr)
    if(ierr/=0) error stop
  elseif (lx3all==1) then
    call hout%read('/I', Incall(:,1), ierr)
    if(ierr/=0) error stop
  else  !< 3d sim, so "/I" is 2D
    call hout%read('/I', Incall, ierr)
    if(ierr/=0) error stop
  endif
  call hout%read('/nullpts', nullptsall, ierr)
  if(ierr/=0) error stop

  call hout%read('/e1', e1all, ierr)
  if(ierr/=0) error stop
  call hout%read('/e2', e2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/e3', e3all, ierr)
  if(ierr/=0) error stop

  call hout%read('/er', erall, ierr)
  if(ierr/=0) error stop
  call hout%read('/etheta', ethetaall, ierr)
  if(ierr/=0) error stop
  call hout%read('/ephi', ephiall, ierr)
  if(ierr/=0) error stop

  call hout%read('/r', rall, ierr)
  if(ierr/=0) error stop
  call hout%read('/theta', thetaall, ierr)
  if(ierr/=0) error stop
  call hout%read('/phi', phiall, ierr)
  if(ierr/=0) error stop

else
  !! 2D with swapped axes

  call hout%read('/x2', x%x3all, ierr)
  if(ierr/=0) error stop
  call hout%read('/x2i', x%x3iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx2b', x%dx3all, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx2h', x%dx3iall, ierr)
  if(ierr/=0) error stop
  !! for a 3D grid this is x2, but now considered x3(all)
  call hout%read('/x3', x%x2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/x3i', x%x2iall, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx3b', x%dx2all, ierr)
  if(ierr/=0) error stop
  call hout%read('/dx3h', x%dx2iall, ierr)
  if(ierr/=0) error stop
  !! formerly x3, now x2

  block
  !> NOTE: workaround for Intel 2020, may not really be a bug
  !> notice this is the only one with negative indices
  !real(wp), dimension(-1:lx1+2,-1:lx3all+2,-1:lx2all+2) :: htmp
  real(wp), allocatable :: htmp(:,:,:)
  allocate(htmp(-1:lx1+2,-1:lx3all+2,-1:lx2all+2))
  !! end workaround

  call hout%read('/h1', htmp, ierr)
  if (ierr/=0) error stop
  x%h1all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
  call hout%read('/h2', htmp, ierr)
  if(ierr/=0) error stop
  !! this would be h3, but with the input structure shape
  x%h3all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
  !! permute the dimensions of the array 3 --> 2, 2 --> 3
  call hout%read('/h3', htmp, ierr)
  if(ierr/=0) error stop
  !! this would be h3, but with the input structure shape
  x%h2all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(1:lx1+1,1:lx3all) :: htmp
  if(lx3all==1) error stop 'lx3 is assumed to be degenerate (?)'
  !! input 2 vs. 3 dimensions swapped from this program
  call hout%read('/h1x1i', htmp, ierr)
  if(ierr/=0) error stop
  x%h1x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/h2x1i', htmp, ierr)
  if(ierr/=0) error stop
  x%h3x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/h3x1i', htmp, ierr)
  if(ierr/=0) error stop
  x%h2x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(1:lx1,1:lx3all+1) :: htmp
  call hout%read('/h1x2i', htmp, ierr)
  if(ierr/=0) error stop
  x%h1x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
  !! Note also that the x2 interface from teh input file is x3i in this simulation
  call hout%read('/h2x2i', htmp, ierr)
  if(ierr/=0) error stop
  x%h3x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
  call hout%read('/h3x2i', htmp, ierr)
  if(ierr/=0) error stop
  x%h2x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
  end block

  block
  real(wp), dimension(1:lx1,1:lx3all,1:lx2all+1) :: htmp
  call hout%read('/h1x3i', htmp, ierr)
  if(ierr/=0) error stop
  x%h1x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  call hout%read('/h2x3i', htmp, ierr)
  if(ierr/=0) error stop
  x%h3x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  call hout%read('/h3x3i', htmp, ierr)
  if(ierr/=0) error stop
  x%h2x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(lx1,lx3all) :: htmp
  call hout%read('/gx1', htmp, ierr)
  if(ierr/=0) error stop
  g1all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/gx2', htmp, ierr)
  if(ierr/=0) error stop
  g3all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/gx3', htmp, ierr)
  if(ierr/=0) error stop
  g2all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])

  call hout%read('/alt', htmp, ierr)
  if(ierr/=0) error stop
  altall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/glat', htmp, ierr)
  if(ierr/=0) error stop
  glatall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/glon', htmp, ierr)
  if(ierr/=0) error stop
  glonall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])

  call hout%read('/Bmag', htmp, ierr)
  if(ierr/=0) error stop
  Bmagall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(lx3all) :: htmp
  call hout%read('/I', htmp, ierr)
  if(ierr/=0) error stop
  Incall = reshape(htmp,[lx2all,lx3all],order=[2,1])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(lx1,lx3all) :: htmp
  call hout%read('/nullpts', htmp, ierr)
  if(ierr/=0) error stop
  nullptsall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  end block

  block
  real(wp), dimension(lx1,lx3all,lx2all,3) :: htmp
  call hout%read('/e1', htmp, ierr)
  if(ierr/=0) error stop
  e1all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  !> swap the x2/x3 unit vectors
  call hout%read('/e2', htmp, ierr)
  if(ierr/=0) error stop
  e3all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  call hout%read('/e3', htmp, ierr)
  if(ierr/=0) error stop
  e2all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  call hout%read('/er', htmp, ierr)
  if(ierr/=0) error stop
  erall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  call hout%read('/etheta', htmp, ierr)
  if(ierr/=0) error stop
  ethetaall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  call hout%read('/ephi', htmp, ierr)
  if(ierr/=0) error stop
  ephiall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  end block

  block
  !! htmp doesn't include the degenerate dimension
  real(wp), dimension(lx1,lx3all) :: htmp
  call hout%read('/r', htmp, ierr)
  if(ierr/=0) error stop
  rall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/theta', htmp, ierr)
  if(ierr/=0) error stop
  thetaall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  call hout%read('/phi', htmp, ierr)
  if(ierr/=0) error stop
  phiall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  end block
endif

call hout%finalize(ierr)
if(ierr/=0) error stop

x%rall = rall
x%thetaall = thetaall
x%phiall = phiall

end procedure get_grid3

end submodule readgrid_hdf5