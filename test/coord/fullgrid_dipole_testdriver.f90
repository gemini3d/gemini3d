program fullgrid_dipole_testdriver

use filesystem, only : mkdir
use phys_consts, only: wp
use meshobj_dipole, only : dipolemesh

implicit none (type, external)

integer, parameter :: lq = 44 + 4, lp = 32 + 4, lphi = 28 + 4
!! +4 for ghost cells
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
! from tohoku20113D_lowres_3Dneu
real(wp), dimension(2), parameter :: qlims=[-0.5340405,0.5340405]
real(wp), dimension(2), parameter :: plims=[1.2509838,1.4372374]
real(wp), dimension(2), parameter :: philims=[3.6126509,3.7240195]
integer :: iq,ip,iphi, i
real(wp), dimension(:,:,:), allocatable :: proj

character(:), allocatable :: path
character(1000) :: argv


! define a grid, in reality this would be pull in from a file
q=[(qlims(1) + (qlims(2)-qlims(1))/(lq-1)*(iq-1),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1))/(lp-1)*(ip-1),ip=1,lp)]
phi=[(philims(1) + (philims(2)-philims(1))/(lphi-1)*(iphi-1),iphi=1,lphi)]

! test min/max coordinate limits
!print*, qlims
!print*, minval(q),maxval(q)

! oddly the destructor does not get called when the program unit terminates; however by
!  putting the variable inside the block we cause it to go out of scope before the program
!  ends and that indeed causes the destructor to get triggered (so we can test it)
!!do while (.true.)
block
type(dipolemesh) :: x

!!!! grid setup and init
! grid spec.
print*, 'fullgrid_testdriver:  Defining curvilinear coordinates...'
call x%set_coords(q,p,phi,p,phi)

! allocations
print*, 'fullgrid_testdriver:  Allocating space for coordinate-specific arrays...'
call x%init()

! call grid generation for this grid def.
print*, 'fullgrid_testdriver:  Calling dipole mesh constructor...'
call x%make()
!!!! end grid setup and init

! check variable allocation and set status
if(.not. (x%xi_alloc_status .and. x%dxi_alloc_status .and. x%difflen_alloc_status .and. x%null_alloc_status .and. &
  x%geog_set_status)) error stop "failed to allocate"

! now do some basic sanity checks
print '(a)', 'fullgrid_testdriver:  Starting basic checks...'
print '(a,1x,i0)', 'fullgrid_testdriver:  grid type...', x%gridflag

print'(a,1x,2F14.3)', 'fullgrid_testdriver, q:', minval(x%q), maxval(x%q)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, p:', minval(x%p), maxval(x%p)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, phi:', minval(x%phi), maxval(x%phi)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, er:', minval(x%er), maxval(x%er)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, etheta:', minval(x%etheta), maxval(x%etheta)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ephi:', minval(x%ephi), maxval(x%ephi)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, eq:', minval(x%eq), maxval(x%eq)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, ep:', minval(x%ep), maxval(x%ep)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, Bmag (nT):', minval(x%Bmag)*1e9, maxval(x%Bmag)*1e9
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gq:', minval(x%gq), maxval(x%gq)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gp:', minval(x%gp), maxval(x%gp)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, gphi:', minval(x%gphi), maxval(x%gphi)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, I:', minval(x%I), maxval(x%I)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, glon:', minval(x%glon), maxval(x%glon)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, glat:', minval(x%glat), maxval(x%glat)
print'(a,1x,2F14.3)', 'fullgrid_testdriver, alt:', minval(x%alt), maxval(x%alt)


allocate(proj(1:lq-4,1:lp-4,1:lphi-4))
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
deallocate(proj)

! test number of null grid points
print*, ' fullgrid_testdriver, number of null grid points:  ',size(x%inull,1)

!> optionally, write out the grid data to a file
if (command_argument_count() >= 1) then
  call get_command_argument(1, argv, status=i)
  if (i /= 0) error stop "could not get user file write path"
  path = trim(argv)
  call mkdir(path)
  print '(a)', ' fullgrid_testdriver, writing grid coords. to: ' // path
  call x%writegrid(path,0)
  call x%writegridall(path,1)
endif

end block
!!end do

end program fullgrid_dipole_testdriver
