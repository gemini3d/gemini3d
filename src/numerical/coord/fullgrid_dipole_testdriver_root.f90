program fullgrid_dipole_testdriver_root

use phys_consts, only: wp
use meshobj_dipole, only : dipolemesh

implicit none

integer, parameter :: lq=384+4,lp=96+4,lphi=64+4
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
real(wp) :: minchkvar,maxchkvar
real(wp), dimension(1:lq-4,1:lp-4,1:lphi-4) :: proj
character(:), allocatable :: path    !use auto-allocation feature
real(wp), dimension(lq-4,2*(lp-4),2*(lphi-4)) :: tmp=0._wp
real(wp), dimension(lq-4+1,2*(lp-4),2*(lphi-4)) :: tmpghost1=0._wp
real(wp), dimension(lq-4,2*(lp-4)+1,2*(lphi-4)) :: tmpghost2=0._wp
real(wp), dimension(lq-4,2*(lp-4),2*(lphi-4)+1) :: tmpghost3=0._wp
real(wp), dimension(-1:(lq-4)+2,-1:2*(lp-4)+2,-1:2*(lphi-4)+2) :: tmpghostall=0._wp


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
!!do while (.true.)
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
                tmp,tmp,tmp, &
                tmp,tmp,tmp)

!!!! end grid setup and init

! check variable allocation and set status
print*, "fullgrid_testdriver:  allocation statuses..."
print*, x%xi_alloc_status,x%dxi_alloc_status,x%dxi_alloc_status_root,x%difflen_alloc_status,x%null_alloc_status,x%geog_set_status
print*, x%coord_alloc_status,x%coord_alloc_status_root

end block
!!end do

end program fullgrid_dipole_testdriver_root
