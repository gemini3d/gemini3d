program fullgrid_cartesian_testdriver

use pathlib, only : mkdir
use phys_consts, only: wp
use meshobj_cart, only : cartmesh

implicit none (type, external)

integer, parameter :: lz=98+4,lx=128+4,ly=144+4     !+4 for ghost cells
real(wp), parameter :: glonctr=207.7_wp,glatctr=65.8_wp
real(wp), dimension(2), parameter :: zlims=[6.2062e+04,9.6794e+05]
real(wp), dimension(2), parameter :: xlims=[-1.5215e+06,1.5215e+06]
real(wp), dimension(2), parameter :: ylims=[-2.0559e+05,2.0559e+05]
real(wp), dimension(lz) :: z
real(wp), dimension(lx) :: xcart
real(wp), dimension(ly) :: y
integer :: iz,ix,iy
real(wp) :: minchkvar,maxchkvar
real(wp), dimension(1:lz-4,1:lx-4,1:ly-4) :: proj
character(:), allocatable :: path                  !use auto-allocation feature


! define a grid, this will include ghost cells
z=[(zlims(1) + (zlims(2)-zlims(1))/(lz-1)*(iz-1),iz=1,lz)]
xcart=[(xlims(1) + (xlims(2)-xlims(1))/(lx-1)*(ix-1),ix=1,lx)]
y=[(ylims(1) + (ylims(2)-ylims(1))/(ly-1)*(iy-1),iy=1,ly)]

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
print*, "fullgrid_testdriver:  allocation statuses..."
print*, x%xi_alloc_status,x%dxi_alloc_status,x%difflen_alloc_status,x%null_alloc_status,x%geog_set_status

! now do some basic sanity checks
print*, 'fullgrid_testdriver:  Starting basic checks...'
print*, 'fullgrid_testdriver:  grid type...',x%gridflag
minchkvar=minval(x%z); maxchkvar=maxval(x%z);
print*, ' fullgrid_testdriver, z:  ',minchkvar,maxchkvar
minchkvar=minval(x%x); maxchkvar=maxval(x%x);
print*, ' fullgrid_testdriver, x:  ',minchkvar,maxchkvar
minchkvar=minval(x%y); maxchkvar=maxval(x%y);
print*, ' fullgrid_testdriver, y:  ',minchkvar,maxchkvar
minchkvar=minval(x%r); maxchkvar=maxval(x%r);
print*, ' fullgrid_testdriver, r:  ',minchkvar,maxchkvar
minchkvar=minval(x%theta); maxchkvar=maxval(x%theta);
print*, ' fullgrid_testdriver, theta:  ',minchkvar,maxchkvar
minchkvar=minval(x%phi); maxchkvar=maxval(x%phi);
print*, ' fullgrid_testdriver, phi:  ',minchkvar,maxchkvar
minchkvar=minval(x%er); maxchkvar=maxval(x%er);
print*, ' fullgrid_testdriver, er:  ',minchkvar,maxchkvar
minchkvar=minval(x%etheta); maxchkvar=maxval(x%etheta);
print*, ' fullgrid_testdriver, etheta:  ',minchkvar,maxchkvar
minchkvar=minval(x%ephi); maxchkvar=maxval(x%ephi);
print*, ' fullgrid_testdriver, ephi:  ',minchkvar,maxchkvar
minchkvar=minval(x%ez); maxchkvar=maxval(x%ez);
print*, ' fullgrid_testdriver, ez:  ',minchkvar,maxchkvar
minchkvar=minval(x%ex); maxchkvar=maxval(x%ex);
print*, ' fullgrid_testdriver, ex:  ',minchkvar,maxchkvar
minchkvar=minval(x%ey); maxchkvar=maxval(x%ey);
print*, ' fullgrid_testdriver, ey:  ',minchkvar,maxchkvar
minchkvar=minval(x%Bmag); maxchkvar=maxval(x%Bmag);
print*, ' fullgrid_testdriver, Bmag (nT):  ',minchkvar*1e9,maxchkvar*1e9
minchkvar=minval(x%gz); maxchkvar=maxval(x%gz);
print*, ' fullgrid_testdriver, gz:  ',minchkvar,maxchkvar
minchkvar=minval(x%gx); maxchkvar=maxval(x%gx);
print*, ' fullgrid_testdriver, gx:  ',minchkvar,maxchkvar
minchkvar=minval(x%gy); maxchkvar=maxval(x%gy);
print*, ' fullgrid_testdriver, gy:  ',minchkvar,maxchkvar
minchkvar=minval(x%I); maxchkvar=maxval(x%I);
print*, ' fullgrid_testdriver, I:  ',minchkvar,maxchkvar
minchkvar=minval(x%glon); maxchkvar=maxval(x%glon)
print*, ' fullgrid_testdriver, glon:  ',minchkvar,maxchkvar
minchkvar=minval(x%glat); maxchkvar=maxval(x%glat)
print*, ' fullgrid_testdriver, glat:  ',minchkvar,maxchkvar
minchkvar=minval(x%alt); maxchkvar=maxval(x%alt)
print*, ' fullgrid_testdriver, alt:  ',minchkvar,maxchkvar

! test number of null grid points
print*, ' fullgrid_testdriver, number of null grid points:  ',size(x%inull,1)

! write out the grid data to a file
path='./cartesian/'
call mkdir(path)
print*, ' fullgrid_testdriver, writing grid coords. to:  ',path
call x%writegrid(path,0)
call x%writegridall(path,1)
end block
!!end do

end program fullgrid_cartesian_testdriver
