program temperature_floor
use phys_consts, only: wp,lsp,kB,gammas
use meshobj_cart, only: cartmesh
use grid, only: set_total_grid_sizes,set_subgrid_sizes
use multifluid, only: rhoe2T,T2rhoe
use transport_audit, only: audit_transport_init,audit_step,audit_mass_start,audit_mass_finish
implicit none(type,external)
type(cartmesh) :: x
real(wp) :: ns(-1:4,-1:4,-1:4,lsp),rhoe(-1:4,-1:4,-1:4,lsp),Ts(-1:4,-1:4,-1:4,lsp)
real(wp) :: expected
integer :: s,unit
character(4096) :: directory
call get_command_argument(1,directory)
call set_total_grid_sizes(2,2,2);call set_subgrid_sizes(2,2)
x%lx1=2;x%lx2=2;x%lx3=2
allocate(x%h1(-1:4,-1:4,-1:4),x%h2(-1:4,-1:4,-1:4),x%h3(-1:4,-1:4,-1:4))
allocate(x%dx1i(2),x%dx2i(2),x%dx3i(2),x%nullpts(2,2,2),x%inull(1,3))
x%h1=1;x%h2=1;x%h3=1;x%dx1i=[1._wp,2._wp];x%dx2i=[3._wp,4._wp];x%dx3i=[5._wp,6._wp]
x%nullpts=.false.;x%nullpts(2,2,2)=.true.;x%lnull=1;x%inull(1,:)=[2,2,2]
ns=1e10_wp;ns(:,:,:,lsp)=6e10_wp
do s=1,lsp
  rhoe(:,:,:,s)=ns(:,:,:,s)*kB*50._wp/(gammas(s)-1)
enddo
call audit_transport_init(trim(directory),0)
call audit_step(0._wp,1._wp)
call audit_mass_start(ns,x)
call rhoe2T(ns,rhoe,Ts)
if(any(Ts/=100._wp)) error stop 'Temperature conversion floor not applied'
call T2rhoe(ns,Ts,rhoe)
do s=1,lsp
  ! Total volume is 3*7*11 minus the 2*4*6 null-cell volume: 183.
  expected=183._wp*ns(1,1,1,s)*kB*50._wp/(gammas(s)-1)
  print '(I0,1X,ES24.16E3)',s,expected
enddo
call audit_mass_finish(ns)
open(newunit=unit,file=trim(directory)//'/temperature.bin',access='stream',form='unformatted',status='replace')
write(unit) Ts,rhoe
close(unit)
end program
