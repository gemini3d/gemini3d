program fullgrid_testdriver

use phys_consts, only: wp
use meshobj_dipole, only : dipolemesh

implicit none

integer, parameter :: lq=256,lp=192,lphi=128
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
real(wp), dimension(2), parameter :: qlims=[-0.5851937,0.5851937]
real(wp), dimension(2), parameter :: plims=[1.2053761,1.5820779]
real(wp), dimension(2), parameter :: philims=[2.0,2.5]
integer :: iq,ip,iphi
real(wp) :: minchkvar,maxchkvar
real(wp), dimension(1:lq-4,1:lp-4,1:lphi-4) :: proj
character(:), allocatable :: path    !use auto-allocation feature


! define a grid, in reality this would be pull in from a file
q=[(qlims(1) + (qlims(2)-qlims(1))/lq*(iq-1),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1))/lp*(ip-1),ip=1,lp)]
phi=[(philims(1) + (philims(2)-philims(1))/lphi*(iphi-1),iphi=1,lphi)]

! oddly the destructor does not get called when the program unit terminates; however by
!  putting the variable inside the block we cause it to go out of scope before the program
!  ends and that indeed causes the destructor to get triggered (so we can test it)
block
type(dipolemesh) :: x


!!!! grid setup and init
! grid spec.
print*, 'fullgrid_testdriver:  Defining curvilinear coordinates...'
call x%set_coords(q,p,phi,p,phi)

! allocations
print*, 'fullgrid_testdriver:  Allocating space for coordinate-specific arrays...'
call x%init_dipolemesh()

! call grid generation for this grid def.
print*, 'fullgrid_testdriver:  Calling dipole mesh constructor...'
call x%make_dipolemesh()

! now some generic methods that can be called once coordinate data are filled in
!print*, 'fullgrid_testdriver:  Setting generic mesh object grid metadata'
!call x%calc_difflengths()
!call x%calc_inull()
!call x%calc_gridflag()
!!!! end grid setup and init

! check variable allocation and set status
print*, "fullgrid_testdriver:  allocation statuses..."
print*, x%xi_alloc_status,x%dxi_alloc_status,x%difflen_alloc_status,x%null_alloc_status,x%geog_set_status

! now do some basic sanity checks
print*, 'fullgrid_testdriver:  Starting basic checks...'
print*, 'fullgrid_testdriver:  grid type...',x%gridflag
minchkvar=minval(x%q); maxchkvar=maxval(x%q);
print*, ' fullgrid_testdriver, q:  ',minchkvar,maxchkvar
minchkvar=minval(x%p); maxchkvar=maxval(x%p);
print*, ' fullgrid_testdriver, p:  ',minchkvar,maxchkvar
minchkvar=minval(x%phi); maxchkvar=maxval(x%phi);
print*, ' fullgrid_testdriver, phi:  ',minchkvar,maxchkvar
minchkvar=minval(x%er); maxchkvar=maxval(x%er);
print*, ' fullgrid_testdriver, er:  ',minchkvar,maxchkvar
minchkvar=minval(x%etheta); maxchkvar=maxval(x%ephi);
print*, ' fullgrid_testdriver, etheta:  ',minchkvar,maxchkvar
minchkvar=minval(x%ephi); maxchkvar=maxval(x%ephi);
print*, ' fullgrid_testdriver, ephi:  ',minchkvar,maxchkvar
minchkvar=minval(x%eq); maxchkvar=maxval(x%eq);
print*, ' fullgrid_testdriver, eq:  ',minchkvar,maxchkvar
minchkvar=minval(x%ep); maxchkvar=maxval(x%ep);
print*, ' fullgrid_testdriver, ep:  ',minchkvar,maxchkvar
minchkvar=minval(x%Bmag); maxchkvar=maxval(x%Bmag);
print*, ' fullgrid_testdriver, Bmag (nT):  ',minchkvar*1e9,maxchkvar*1e9
minchkvar=minval(x%gq); maxchkvar=maxval(x%gq);
print*, ' fullgrid_testdriver, gq:  ',minchkvar,maxchkvar
minchkvar=minval(x%gp); maxchkvar=maxval(x%gp);
print*, ' fullgrid_testdriver, gp:  ',minchkvar,maxchkvar
minchkvar=minval(x%gphi); maxchkvar=maxval(x%gphi);
print*, ' fullgrid_testdriver, gphi:  ',minchkvar,maxchkvar
minchkvar=minval(x%I); maxchkvar=maxval(x%I);
print*, ' fullgrid_testdriver, I:  ',minchkvar,maxchkvar
minchkvar=minval(x%glon); maxchkvar=maxval(x%glon)
print*, ' fullgrid_testdriver, glon:  ',minchkvar,maxchkvar
minchkvar=minval(x%glat); maxchkvar=maxval(x%glat)
print*, ' fullgrid_testdriver, glat:  ',minchkvar,maxchkvar
minchkvar=minval(x%alt); maxchkvar=maxval(x%alt)
print*, ' fullgrid_testdriver, alt:  ',minchkvar,maxchkvar


! check orthogonality of the basis vectors
proj=sum(x%eq*x%ep,dim=4)
if (any(abs(proj)>1e-4)) error stop '  eq,ep not ortho!!!'
proj=sum(x%eq*x%ephi,dim=4)
if (any(abs(proj)>1e-4)) error stop '  eq,ephi not ortho!!!'
proj=sum(x%ep*x%ephi,dim=4)
if (any(abs(proj)>1e-4)) error stop '  ep,ephi not ortho!!!'
proj=sum(x%er*x%etheta,dim=4)
if (any(proj>1e-4)) error stop '  er,etheta not ortho!!!'
proj=sum(x%er*x%ephi,dim=4)
if (any(proj>1e-4)) error stop '  er,ephi not ortho!!!'
proj=sum(x%etheta*x%ephi,dim=4)
if (any(proj>1e-4)) error stop '  etheta,ephi not ortho!!!'

! test number of null grid points
print*, ' fullgrid_testdriver, number of null grid points:  ',size(x%inull,1)

! write out the grid data to a file
path='/Users/zettergm/Downloads/'
print*, ' fullgrid_testdriver, writing grid coords. to file...'
call x%writegrid(path,0)
call x%writegridall(path,1)
end block

! deallocate the grid before ending the program
print*, 'fullgrid_testdriver:  exiting program; destructor should have been called...'

end program fullgrid_testdriver
