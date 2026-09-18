program continuity_floor
use phys_consts, only: wp,lsp,mindens
use meshobj_cart, only: cartmesh
use grid, only: set_total_grid_sizes,set_subgrid_sizes
use multifluid, only: clean_param
use transport_audit, only: audit_transport_init,audit_step,audit_mass_start,audit_mass_finish
implicit none(type,external)
type(cartmesh) :: x
real(wp) :: ns(-1:4,-1:4,-1:4,lsp)
character(4096) :: directory
call get_command_argument(1,directory)
call set_total_grid_sizes(2,2,2);call set_subgrid_sizes(2,2)
x%lx1=2;x%lx2=2;x%lx3=2
allocate(x%h1(-1:4,-1:4,-1:4),x%h2(-1:4,-1:4,-1:4),x%h3(-1:4,-1:4,-1:4))
allocate(x%dx1i(2),x%dx2i(2),x%dx3i(2),x%nullpts(2,2,2),x%inull(1,3))
x%h1=1;x%h2=1;x%h3=1;x%dx1i=1;x%dx2i=1;x%dx3i=1
x%nullpts=.false.;x%nullpts(2,2,2)=.true.;x%lnull=1;x%inull(1,:)=[2,2,2]
ns=.5_wp;ns(:,:,:,7)=3;mindens=2
call audit_transport_init(trim(directory),0)
call audit_step(0._wp,1._wp)
call audit_mass_start(ns,x)
call clean_param(x,1,ns)
call audit_mass_finish(ns)
end program
