program audit_density_floor
use phys_consts, only: wp,lsp,mindens,mindensnull,mindensdiv
use multifluid, only: clean_param
use meshobj_cart, only: cartmesh
use grid, only: set_total_grid_sizes,set_subgrid_sizes
implicit none
real(wp)::ns(-1:4,-1:3,-1:3,lsp),floor_value
real(wp),parameter::floors(3)=[0.01_wp,1000._wp,5000._wp]
type(cartmesh)::x
integer::i
call set_total_grid_sizes(2,1,1)
call set_subgrid_sizes(1,1)
x%lnull=1
allocate(x%inull(1,3));x%inull(1,:)=[2,1,1]
do i=1,size(floors)
  floor_value=floors(i);mindens=floor_value
  ns=0
  call clean_param(x,1,ns)
  if(any(ns(1,1,1,1:6)/=floor_value))error stop 'configured ion floor was overridden'
  if(abs(ns(1,1,1,7)-6*floor_value)>1e-12_wp)error stop 'physical-cell charge neutrality'
  if(ns(2,1,1,1)/=mindensnull*1e-2_wp)error stop 'null-cell fill changed'
  if(any(ns(-1,:,:,:)/=mindensdiv))error stop 'ghost fill changed'
enddo
deallocate(x%inull)
print *, 'Configured density floor and ESF opt-in passed'
end program
