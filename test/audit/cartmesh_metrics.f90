module audit_checked_cartmesh
use phys_consts, only: wp,Re
use meshobj_cart, only: cartmesh
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
type, extends(cartmesh) :: checked_cartmesh
contains
  procedure, nopass :: calc_h1=>checked_metric
  procedure, nopass :: calc_h2=>checked_metric
  procedure, nopass :: calc_h3=>checked_metric
end type
contains
function checked_metric(coord1,coord2,coord3) result(hval)
  real(wp), pointer, intent(in) :: coord1(:,:,:),coord2(:,:,:),coord3(:,:,:)
  real(wp) :: hval(lbound(coord1,1):ubound(coord1,1),lbound(coord1,2):ubound(coord1,2), &
    lbound(coord1,3):ubound(coord1,3))
  if (.not.associated(coord1) .or. .not.associated(coord2) .or. .not.associated(coord3)) &
    error stop 'Cartesian metric received deallocated coordinates'
  if (any(lbound(coord1)/=-1)) error stop 'Cartesian coordinate ghost bounds lost'
  if (.not.all(ieee_is_finite(coord1)) .or. .not.all(ieee_is_finite(coord2)) .or. &
      .not.all(ieee_is_finite(coord3))) error stop 'Nonfinite Cartesian coordinates'
  if (any(coord1<Re+70000._wp)) error stop 'Cartesian radial coordinates corrupted'
  hval=1
end function
end module

program audit_cartmesh_metrics
use audit_checked_cartmesh, only: checked_cartmesh
use phys_consts, only: wp
implicit none
integer :: iteration
do iteration=1,4
  call check_mesh(5,4,3)
  call check_mesh(5,1,4)
  call check_mesh(5,4,1)
enddo
print *, 'Cartesian metrics retain live coordinates, ghost bounds and singleton dimensions'
contains
subroutine check_mesh(n1,n2,n3)
  integer, intent(in) :: n1,n2,n3
  type(checked_cartmesh) :: x
  real(wp) :: z(-1:n1+2),east(-1:n2+2),north(-1:n3+2)
  integer :: i
  z=[(100000._wp+10000._wp*i,i=-1,n1+2)]
  east=[(1000._wp*i,i=-1,n2+2)]
  north=[(2000._wp*i,i=-1,n3+2)]
  call x%set_coords(z,east,north,east,north)
  call x%set_center(0._wp,45._wp)
  call x%init()
  call x%make()
  if (any(shape(x%h1)/=[n1+4,n2+4,n3+4])) error stop 'Cartesian metric shape'
  if (any(x%h1/=1) .or. any(x%h2/=1) .or. any(x%h3/=1)) error stop 'Cartesian center metrics'
  if (any(x%h1x1i/=1) .or. any(x%h2x1i/=1) .or. any(x%h3x1i/=1)) error stop 'Cartesian x1 interfaces'
  if (any(x%h1x2i/=1) .or. any(x%h2x2i/=1) .or. any(x%h3x2i/=1)) error stop 'Cartesian x2 interfaces'
  if (any(x%h1x3i/=1) .or. any(x%h2x3i/=1) .or. any(x%h3x3i/=1)) error stop 'Cartesian x3 interfaces'
end subroutine
end program
