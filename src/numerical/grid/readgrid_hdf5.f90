submodule (grid) readgrid_hdf5

use phys_consts, only: debug
use h5fortran, only: hdf5_file

implicit none (type, external)

contains


module procedure get_grid3_coords_hdf5
  type(hdf5_file) :: hf
  character(:), allocatable :: fn

  if (index(path, 'simgrid.h5') /= 0) then
    fn = path
  else
    fn = path // '/simgrid.h5'
  endif
  if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid:', fn

  call hf%open(fn, action='r')

  !> reads common to 2D and 3D

  !> 1-D variables
  call hf%read('/x1', x1)
  call hf%read('/x2', x2all)
  call hf%read('/x3', x3all)

  if (.not. (hf%exist('/glonctr') .and. &
             hf%exist('/glatctr'))) error stop "please specify /glonctr and /glatctr in " // fn

  call hf%read('/glonctr',glonctr)
  call hf%read('/glatctr',glatctr)

  call hf%close()
end procedure get_grid3_coords_hdf5


!module procedure get_grid3_hdf5
!
!  type(hdf5_file) :: hout
!  character(:), allocatable :: fn
!
!  if (index(path, 'simgrid.h5') /= 0) then
!    fn = path
!  else
!    fn = path // '/simgrid.h5'
!  endif
!  if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid:', fn
!
!  call hout%open(fn, action='r')
!
!  !> reads common to 2D and 3D
!
!  !> 1-D variables
!  call hout%read('/x1', x%x1)
!  call hout%read('/x1i', x%x1i)
!  call hout%read('/dx1b', x%dx1)
!  call hout%read('/dx1h', x%dx1i)
!
!  if (flagswap/=1) then
!    !! normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
!
!    call hout%read('/x2', x%x2all)
!    call hout%read('/x2i', x%x2iall)
!    call hout%read('/dx2b', x%dx2all)
!    call hout%read('/dx2h', x%dx2iall)
!
!    call hout%read('/x3', x%x3all)
!    call hout%read('/x3i', x%x3iall)
!    call hout%read('/dx3b', x%dx3all)
!    call hout%read('/dx3h', x%dx3iall)
!
!    !> 3D variables
!    call hout%read('/h1', x%h1all)
!    call hout%read('/h2', x%h2all)
!    call hout%read('/h3', x%h3all)
!
!    call hout%read('/h1x1i', x%h1x1iall)
!    call hout%read('/h2x1i', x%h2x1iall)
!    call hout%read('/h3x1i', x%h3x1iall)
!
!    call hout%read('/h1x2i', x%h1x2iall)
!    call hout%read('/h2x2i', x%h2x2iall)
!    call hout%read('/h3x2i', x%h3x2iall)
!
!    call hout%read('/h1x3i', x%h1x3iall)
!    call hout%read('/h2x3i', x%h2x3iall)
!    call hout%read('/h3x3i', x%h3x3iall)
!
!    call hout%read('/gx1', g1all)
!    call hout%read('/gx2', g2all)
!    call hout%read('/gx3', g3all)
!
!    call hout%read('/alt', x%altall)
!    call hout%read('/glat', glatall)
!    call hout%read('/glon', x%glonall)
!
!    call hout%read('/Bmag', x%Bmagall)
!    call hout%read('/I', Incall)
!    call hout%read('/nullpts', nullptsall)
!
!    call hout%read('/e1', e1all)
!    call hout%read('/e2', e2all)
!    call hout%read('/e3', e3all)
!
!    call hout%read('/er', erall)
!    call hout%read('/etheta', ethetaall)
!    call hout%read('/ephi', ephiall)
!
!    call hout%read('/r', rall)
!    call hout%read('/theta', thetaall)
!    call hout%read('/phi', phiall)
!
!  else
!    !! 2D with swapped axes
!
!    call hout%read('/x2', x%x3all)
!    call hout%read('/x2i', x%x3iall)
!    call hout%read('/dx2b', x%dx3all)
!    call hout%read('/dx2h', x%dx3iall)
!    !! for a 3D grid this is x2, but now considered x3(all)
!    call hout%read('/x3', x%x2all)
!    call hout%read('/x3i', x%x2iall)
!    call hout%read('/dx3b', x%dx2all)
!    call hout%read('/dx3h', x%dx2iall)
!    !! formerly x3, now x2
!
!    block
!    !> NOTE: workaround for Intel 2020, may not really be a bug
!    !> notice this is the only one with negative indices
!    !real(wp), dimension(-1:lx1+2,-1:lx3all+2,-1:lx2all+2) :: htmp
!    real(wp), allocatable :: htmp(:,:,:)
!    allocate(htmp(-1:lx1+2,-1:lx3all+2,-1:lx2all+2))
!    !! end workaround
!
!    call hout%read('/h1', htmp)
!    x%h1all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
!    call hout%read('/h2', htmp)
!    !! this would be h3, but with the input structure shape
!    x%h3all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
!    !! permute the dimensions of the array 3 --> 2, 2 --> 3
!    call hout%read('/h3', htmp)
!    !! this would be h3, but with the input structure shape
!    x%h2all = reshape(htmp, [lx1+4,lx2all+4,lx3all+4], order=[1,3,2])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(1:lx1+1,1:lx3all, 1) :: htmp
!    if(lx3all==1) error stop 'lx3 is assumed to be degenerate'
!    !! input 2 vs. 3 dimensions swapped from this program
!    call hout%read('/h1x1i', htmp)
!    x%h1x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/h2x1i', htmp)
!    x%h3x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/h3x1i', htmp)
!    x%h2x1iall = reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(1:lx1,1:lx3all+1, 1) :: htmp
!    call hout%read('/h1x2i', htmp)
!    x%h1x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
!    !! Note also that the x2 interface from teh input file is x3i in this simulation
!    call hout%read('/h2x2i', htmp)
!    x%h3x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
!    call hout%read('/h3x2i', htmp)
!    x%h2x3iall = reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(1:lx1,1:lx3all,1:lx2all+1) :: htmp
!    call hout%read('/h1x3i', htmp)
!    x%h1x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    call hout%read('/h2x3i', htmp)
!    x%h3x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    call hout%read('/h3x3i', htmp)
!    x%h2x2iall = reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(lx1, lx3all, 1) :: htmp
!    call hout%read('/gx1', htmp)
!    g1all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/gx2', htmp)
!    g3all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/gx3', htmp)
!    g2all = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!
!    call hout%read('/alt', htmp)
!    x%altall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/glat', htmp)
!    glatall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/glon', htmp)
!    x%glonall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!
!    call hout%read('/Bmag', htmp)
!    x%Bmagall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(lx3all, 1) :: htmp
!    call hout%read('/I', htmp)
!    Incall = reshape(htmp,[lx2all,lx3all],order=[2,1])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(lx1, lx3all, 1) :: htmp
!    call hout%read('/nullpts', htmp)
!    nullptsall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(lx1,lx3all,lx2all,3) :: htmp
!    call hout%read('/e1', htmp)
!    e1all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!
!    !> swap the x2/x3 unit vectors
!    call hout%read('/e2', htmp)
!    e3all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    call hout%read('/e3', htmp)
!    e2all = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    call hout%read('/er', htmp)
!    erall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    call hout%read('/etheta', htmp)
!    ethetaall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    call hout%read('/ephi', htmp)
!    ephiall = reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    end block
!
!    block
!    !! htmp includes degenerate dimension
!    real(wp), dimension(lx1, lx3all, 1) :: htmp
!    call hout%read('/r', htmp)
!    rall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/theta', htmp)
!    thetaall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    call hout%read('/phi', htmp)
!    phiall = reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!  endif
!
!  call hout%close()
!
!  x%rall = rall
!  x%thetaall = thetaall
!  x%phiall = phiall
!
!end procedure get_grid3_hdf5
!
end submodule readgrid_hdf5
