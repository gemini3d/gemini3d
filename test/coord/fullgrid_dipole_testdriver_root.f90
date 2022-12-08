program fullgrid_dipole_testdriver_root

use, intrinsic :: iso_fortran_env, only : int64
use phys_consts, only: wp
use meshobj_dipole, only : dipolemesh

implicit none (type, external)

integer, parameter :: lq = 44 + 4, lp = 32 + 4, lphi = 28 + 4
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(lphi) :: phi
real(wp), dimension(2*(lp-4)+4) :: pall
real(wp), dimension(2*(lphi-4)+4) :: phiall
! from tohoku20113D_lowres_3Dneu
real(wp), dimension(2), parameter :: qlims=[-0.5340405,0.5340405]
real(wp), dimension(2), parameter :: plims=[1.2509838,1.4372374]
real(wp), dimension(2), parameter :: philims=[3.6126509,3.7240195]
integer :: iq,ip,iphi
integer :: ierr, i, N
integer(int64) :: mem_bytes, Bel

character(10) :: argv

real(wp), allocatable, dimension(:,:,:) :: tmp, tmpghost1, tmpghost2, tmpghost3, tmpghostall


allocate(tmp(lq-4,2*(lp-4),2*(lphi-4)), &
tmpghost1(lq-4+1,2*(lp-4),2*(lphi-4)), &
tmpghost2(lq-4,2*(lp-4)+1,2*(lphi-4)), &
tmpghost3(lq-4,2*(lp-4),2*(lphi-4)+1), &
tmpghostall(-1:(lq-4)+2,-1:2*(lp-4)+2,-1:2*(lphi-4)+2), &
stat=ierr)

if(ierr /= 0) error stop "failed to allocate memory. May need to try a smaller test grid"

!> estimate memory used
Bel = storage_size(q, kind=int64) / 8

mem_bytes = Bel * 2 * (6 * size(tmp, kind=int64) + 3 * size(tmpghostall, kind=int64) +  &
  3 * size(tmpghost1, kind=int64) + 3 * size(tmpghost2, kind=int64) + 3 * size(tmpghost3, kind=int64))
print *, "estimated memory used (Megabytes): ", mem_bytes / 1000000

tmp = 0
tmpghost1 = 0
tmpghost2 = 0
tmpghost3 = 0
tmpghostall = 0

N = 1
if(command_argument_count() >= 1) then
  call get_command_argument(1, argv, status=i)
  if (i/=0) error stop "first argument is number of loops"
  read(argv,'(I4)') N
endif

! define a grid, in reality this would be pulled in from a file
q=[(qlims(1) + (qlims(2)-qlims(1))/(lq-1)*(iq-1),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1))/(lp-1)*(ip-1),ip=1,lp)]
phi=[(philims(1) + (philims(2)-philims(1))/(lphi-1)*(iphi-1),iphi=1,lphi)]

! define a global grid, doesn't really matter how this overlaps the subgrid...
pall=[(plims(1) + (plims(2)-plims(1))/(lp-1)*(ip-1),ip=1,2*(lp-4)+4)]
phiall=[(philims(1) + (philims(2)-philims(1))/(lphi-1)*(iphi-1),iphi=1,2*(lphi-4)+4)]

! oddly the destructor does not get called when the program unit terminates; however by
!  putting the variable inside the block we cause it to go out of scope before the program
!  ends and that indeed causes the destructor to get triggered (so we can test it)
main : do i = 1,N
!! we loop N times to help ensure memory isn't leaking,
!! that we free appropriate variables upon destruction

block
type(dipolemesh) :: x

!!!! grid setup and init
! grid spec.
print*, 'fullgrid_testdriver:  Defining curvilinear coordinates...'
call x%set_coords(q,p,phi,pall,phiall)

! allocations
print*, 'fullgrid_testdriver:  Allocating space for coordinate-specific arrays...'
call x%init()

! allocate space for root grid
print*, 'fullgrid_testdriver:  Allocating space for root arrays...'
call x%init_storage_root()

! call grid generation for this grid def.
print*, 'fullgrid_testdriver:  Calling dipole mesh constructor...'
call x%make()

! fill root array data
print*, 'fullgrid_testdriver:  setting root grid data:  ',x%lx1,x%lx2,x%lx3,x%lx2all,x%lx3all
!print*, shape(tmp),shape(tmpghost1),shape(tmpghost2),shape(tmpghost3),shape(tmpghostall)
call x%set_root(tmpghostall,tmpghostall,tmpghostall, &
                tmpghost1,tmpghost1,tmpghost1, &
                tmpghost2,tmpghost2,tmpghost2, &
                tmpghost3,tmpghost3,tmpghost3, &
                tmpghostall,tmpghostall,tmpghostall, &
                tmpghostall,tmp,tmpghostall)
call x%calc_coord_diffs_root()
!!!! end grid setup and init

! check variable allocation and set status
if (.not. x%xi_alloc_status) error stop "xi alloc false"
if (.not. x%dxi_alloc_status) error stop "dxi alloc false"
if (.not. x%dxi_alloc_status_root) error stop "dxi_root alloc false"
if (.not. x%difflen_alloc_status) error stop "difflen alloc false"
if (.not. x%null_alloc_status) error stop "null alloc false"
if (.not. x%geog_set_status) error stop "geog alloc false"
if (.not. x%coord_alloc_status) error stop "coord alloc false"
if (.not. x%coord_alloc_status_root) error stop "coord_root alloc false"

end block

end do main

end program fullgrid_dipole_testdriver_root
