program fullgrid_cartesian_testdriver

use filesystem, only : mkdir
use phys_consts, only: wp
use meshobj_cart, only : cartmesh

implicit none (type, external)

integer, parameter :: lz = 44 + 4, lx = 32 + 4, ly = 28 + 4
!! +4 for ghost cells
real(wp), parameter :: glonctr = 207.7, glatctr = 65.8
real(wp), dimension(2), parameter :: zlims=[6.2062e+04, 9.6794e+05]
real(wp), dimension(2), parameter :: xlims=[-1.5215e+06, 1.5215e+06]
real(wp), dimension(2), parameter :: ylims=[-2.0559e+05, 2.0559e+05]
real(wp), dimension(lz) :: z
real(wp), dimension(lx) :: xcart
real(wp), dimension(ly) :: y
integer :: iz,ix,iy, i

real(wp), dimension(:,:,:), allocatable :: proj

character(:), allocatable :: path
character(1000) :: argv

allocate(proj(1:lz-4,1:lx-4,1:ly-4))

! define a grid, this will include ghost cells
z=[(zlims(1) + (zlims(2)-zlims(1)) / (lz-1)*(iz-1),iz=1,lz)]
xcart=[(xlims(1) + (xlims(2)-xlims(1)) / (lx-1)*(ix-1),ix=1,lx)]
y=[(ylims(1) + (ylims(2)-ylims(1)) / (ly-1)*(iy-1),iy=1,ly)]

! oddly the destructor does not get called when the program unit terminates; however by
!  putting the variable inside the block we cause it to go out of scope before the program
!  ends and that indeed causes the destructor to get triggered (so we can test it)
!!do while (.true.)
block
type(cartmesh) :: x

!!!! grid setup and init
! grid spec.
print*, 'fullgrid_testdriver:  Defining curvilinear coordinates...'
call x%set_coords(z,xcart,y,xcart,y)

print*, 'fullgrid_testdriver:  Setting geographic center of grid...'
call x%set_center(glonctr,glatctr)

! allocations
print*, 'fullgrid_testdriver:  Allocating space for coordinate-specific arrays...'
call x%init()

! call grid generation for this grid def.
print*, 'fullgrid_testdriver:  Calling dipole mesh constructor...'
call x%make()
!!!! end grid setup and init

! check variable allocation and set status
if(.not. all([x%xi_alloc_status,x%dxi_alloc_status,x%difflen_alloc_status,x%null_alloc_status,x%geog_set_status])) &
  error stop "grid allocation failure"

! now do some basic sanity checks
print '(a)', 'fullgrid_testdriver:  Starting basic checks...'
print '(a,1x,i0)', 'fullgrid_testdriver:  grid type...', x%gridflag
print '(a,1x,2F14.3)', 'fullgrid_testdriver, z:',  minval(x%z),maxval(x%z)
print '(a,1x,2F14.3)', 'fullgrid_testdriver, x:',  minval(x%x), maxval(x%x)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, y:',  minval(x%y), maxval(x%y)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, r:',  minval(x%r), maxval(x%r)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, theta:',  minval(x%theta), maxval(x%theta)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, phi:',  minval(x%phi), maxval(x%phi)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, er:',  minval(x%er), maxval(x%er)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, etheta:',  minval(x%etheta), maxval(x%etheta)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ephi:',  minval(x%ephi), maxval(x%ephi)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ez:',  minval(x%ez), maxval(x%ez)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ex:',  minval(x%ex), maxval(x%ex)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ey:',  minval(x%ey), maxval(x%ey)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, Bmag (nT):', minval(x%Bmag)*1e9, maxval(x%Bmag)*1e9
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gz:', minval(x%gz), maxval(x%gz)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gx:', minval(x%gx), maxval(x%gx)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gy:', minval(x%gy), maxval(x%gy)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, I:', minval(x%I), maxval(x%I)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, glon:', minval(x%glon), maxval(x%glon)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, glat:', minval(x%glat), maxval(x%glat)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, alt:', minval(x%alt), maxval(x%alt)

! test number of null grid points
print '(a,1x,i0)', 'fullgrid_testdriver, number of null grid points:', size(x%inull,1)

! write out the grid data to a file
if (command_argument_count() >= 1) then
  call get_command_argument(1, argv, status=i)
  if (i /= 0) error stop "could not get user file write path"
  path = trim(argv)
  call mkdir(path)
  print '(a)', 'fullgrid_testdriver, writing grid coords. to: ' // path
  call x%writegrid(path,0)
  call x%writegridall(path,1)
endif

end block
!!end do

end program fullgrid_cartesian_testdriver
