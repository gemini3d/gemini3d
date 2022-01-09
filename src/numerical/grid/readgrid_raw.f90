submodule (grid) readgrid_raw

use phys_consts, only: debug

implicit none (type, external)

contains


module procedure get_grid3_coords_raw

  integer :: u
  character(:), allocatable :: fn
  integer :: stat

  if (index(path, 'simgrid.dat') /= 0) then
    fn = path
  else
    fn = path // '/simgrid.dat'
  endif
  if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid:', fn

  open(newunit=u, file=fn, status='old', form='unformatted', access='stream', action='read')

  !! normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
  read(u) x1
  read(u) x2all
  read(u) x3all

  !! check whether center location is include in file
  read(u,iostat=stat) glonctr
  if (stat/=0) error stop "please specify /glonctr and /glatctr in " // fn
  read(u) glatctr

  close(u)

end procedure get_grid3_coords_raw


!module procedure get_grid3_raw
!
!  integer :: u
!  character(:), allocatable :: fn
!
!  if (index(path, 'simgrid.dat') /= 0) then
!    fn = path
!  else
!    fn = path // '/simgrid.dat'
!  endif
!  if (debug) print '(A,/,A)', 'READ 3D (B-parallel, B-perp, B-perp) grid:', fn
!
!  open(newunit=u, file=fn, status='old', form='unformatted', access='stream', action='read')
!
!  if (flagswap/=1) then
!    !! normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
!    read(u) x%x1,x%x1i,x%dx1,x%dx1i
!    read(u) x%x2all,x%x2iall,x%dx2all,x%dx2iall
!    read(u) x%x3all,x%x3iall,x%dx3all,x%dx3iall
!    read(u) x%h1all,x%h2all,x%h3all
!    read(u) x%h1x1iall,x%h2x1iall,x%h3x1iall
!    read(u) x%h1x2iall,x%h2x2iall,x%h3x2iall
!    read(u) x%h1x3iall,x%h2x3iall,x%h3x3iall
!
!    read(u) g1all,g2all,g3all
!    read(u) x%altall
!    read(u) glatall,x%glonall
!    read(u) x%Bmagall
!    read(u) Incall
!    read(u) nullptsall
!    read(u) e1all
!    read(u) e2all
!    read(u) e3all
!    read(u) erall
!    read(u) ethetaall
!    read(u) ephiall
!
!    read(u) rall
!    read(u) thetaall
!    read(u) phiall
!  else
!    read(u) x%x1,x%x1i,x%dx1,x%dx1i                !< x1 untouched
!    read(u) x%x3all,x%x3iall,x%dx3all,x%dx3iall    !< for a 3D grid this is x2, but now considered x3(all)
!    read(u) x%x2all,x%x2iall,x%dx2all,x%dx2iall    !< formerly x3, now x2
!
!    block
!    real(wp), dimension(-1:lx1+2,-1:lx3all+2,-1:lx2all+2) :: htmp
!    !! this stores the input metric factors which are swapped x2/x3 vs. what this simulation will use
!    read(u) htmp
!    x%h1all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])
!    read(u) htmp        !< this would be h3, but with the input structure shape
!    x%h3all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])   !< permute the dimensions of the array 3 --> 2, 2 --> 3
!    read(u) htmp        !< this would be h3, but with the input structure shape
!    x%h2all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(1:lx1+1,1:lx3all,1:lx2all) :: htmp
!    !! input 2 vs. 3 dimensions swapped from this program
!    read(u) htmp
!    x%h1x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    x%h3x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    x%h2x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(1:lx1,1:lx3all+1,1:lx2all) :: htmp
!    read(u) htmp
!    x%h1x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])    !Note also that the x2 interface from teh input file is x3i in this simulation
!    read(u) htmp
!    x%h3x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
!    read(u) htmp
!    x%h2x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(1:lx1,1:lx3all,1:lx2all+1) :: htmp
!    read(u) htmp
!    x%h1x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    read(u) htmp
!    x%h3x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    read(u) htmp
!    x%h2x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(lx1,lx3all,lx2all) :: htmp
!    read(u) htmp
!    g1all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    g3all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    g2all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!
!    read(u) htmp
!    x%altall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    glatall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    x%glonall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!
!    read(u) htmp
!    x%Bmagall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(lx3all,lx2all) :: htmp
!    read(u) htmp
!    Incall=reshape(htmp,[lx2all,lx3all],order=[2,1])
!    end block
!
!      !inquire(u, pos=itell)
!    !print *,'file pos before read',itell
!    !read(u) nullptsall
!    !print *,'lx1',lx1,'lx2',lx2,'lx3all',lx3all
!    !print *,'shape(nullptsall)',shape(nullptsall)
!    !inquire(u, pos=itell)
!    !print *,'file pos after read',itell
!  ! FIXME BROKEN!
!    !allocate(nullptsall(lx1,lx2,lx3all))
!    !print *,shape(nullptsall)
!    !stop
!
!    ! FIXME would be like this, but this doesn't work.
!    !allocate(nullptsall(lx1,lx2,lx3all))
!    !read(u) nullptsall
!
!    block
!    real(wp), dimension(lx1,lx3all,lx2all) :: htmp
!    read(u) htmp
!    nullptsall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!
!    block
!    real(wp), dimension(lx1,lx3all,lx2all,3) :: htmp
!    read(u) htmp
!    e1all=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!
!    !! swap the x2/x3 unit vectors
!    read(u) htmp
!    e3all=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    read(u) htmp
!    e2all=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    read(u) htmp
!    erall=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    read(u) htmp
!    ethetaall=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    read(u) htmp
!    ephiall=reshape(htmp,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
!    end block
!
!    block
!    real(wp), dimension(lx1,lx3all,lx2all) :: htmp
!    read(u) htmp
!    rall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    thetaall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    read(u) htmp
!    phiall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
!    end block
!  endif
!
!  close(u)
!
!  x%rall=rall
!  x%thetaall=thetaall
!  x%phiall=phiall
!
!end procedure get_grid3_raw
end submodule readgrid_raw
