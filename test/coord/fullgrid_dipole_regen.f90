program fullgrid_dipole_testdriver

use filesystem, only : mkdir
use phys_consts, only: wp
use grid, only: generate_worker_grid, ungenerate_worker_grid
use meshobj_dipole, only : dipolemesh

implicit none (type, external)

integer, parameter :: lq=384+4,lp=96+4,lphi=64+4
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
! from tohoku20113D_lowres_3Dneu
real(wp), dimension(2), parameter :: qlims=[-0.5340405,0.5340405]
real(wp), dimension(2), parameter :: plims=[1.2509838,1.4372374]
real(wp), dimension(2), parameter :: philims=[3.6126509,3.7240195]
integer :: iq,ip,iphi, i
real(wp) :: minchkvar,maxchkvar
real(wp), dimension(:,:,:), allocatable :: proj

character(:), allocatable :: path
character(1000) :: argv

allocate(proj(1:lq-4,1:lp-4,1:lphi-4))

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
integer :: irepeat

do irepeat=1,1000
  call generate_worker_grid(q,p,phi,p,phi,0._wp,0._wp,x)

  ! check variable allocation and set status
  if(.not. (x%xi_alloc_status .and. x%dxi_alloc_status .and. x%difflen_alloc_status .and. x%null_alloc_status .and. &
    x%geog_set_status)) error stop "failed to allocate"
  
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

  ! deallocate and try again
  call ungenerate_worker_grid(x)
end do


!> optionally, write out the grid data to a file
if (command_argument_count() >= 1) then
  call get_command_argument(1, argv, status=i)
  if (i /= 0) error stop "could not get user file write path"
  path = trim(argv)
  call mkdir(path)
  print*, ' fullgrid_testdriver, writing grid coords. to:  ',path
  call x%writegrid(path,0)
  call x%writegridall(path,1)
endif

end block
!!end do

end program fullgrid_dipole_testdriver
