program physics_native
use phys_consts, only: wp,lsp,ln,kb,ms,gammas
use collisions, only: coulomb_colln,thermal_conduct
use sources, only: srcsEnergy
use calculus, only: ETD_uncoupled
use meshobj_cart, only: cartmesh
use grid, only: set_total_grid_sizes,set_subgrid_sizes,set_gridflag
use diffusion, only: diffusion_prep,backEuler3D,TRBDF23D
use nonlinear_diffusion, only: midpoint_TRBDF23D
use diffusion_features, only: midpoint_conduction_enabled
use multifluid, only: energy_diffusion
implicit none(type,external)
character(4096) :: directory,trace_name
integer :: trace_status
call get_command_argument(1,directory)
if(len_trim(directory)==0) error stop 'output directory required'
call nonlinear_conduction()
call kirchhoff_flux()
call get_environment_variable('GEMINI_PHYSICAL_TR_TRACE',trace_name,status=trace_status)
if(trace_status/=0) call wrapper_checks()
contains

subroutine nonlinear_conduction()
integer,parameter :: n=49
type(cartmesh) :: x
real(wp) :: ns(-1:n+2,-1:3,-1:3),temp(-1:n+2,-1:3,-1:3),old(-1:n+2,-1:3,-1:3),raw(-1:n+2,-1:3,-1:3)
real(wp) :: nn(n,1,1,ln),j1(-1:n+2,-1:3,-1:3),tn(n,1,1)
real(wp) :: lam(n,1,1),beta(n,1,1),a(n,1,1),b(n,1,1),c(n,1,1),d(n,1,1),e(n,1,1),inc(n,1,1,5)
real(wp) :: density,dx,total,dt,cv
integer :: id,window,level,method,i,it,nt,u
character(64) :: name
call set_total_grid_sizes(n,1,1);call set_subgrid_sizes(1,1);call set_gridflag(2)
x%lx1=n;x%lx2=1;x%lx3=1
allocate(x%h1(-1:n+2,-1:3,-1:3),x%h2(-1:n+2,-1:3,-1:3),x%h3(-1:n+2,-1:3,-1:3))
allocate(x%dx1(0:n+2),x%dx1i(n))
x%h1=1;x%h2=1;x%h3=1;dx=1e5_wp/(n-1);x%dx1=dx;x%dx1i=dx
nn=0;j1=0;tn=600
do id=1,2
  density=1e5_wp
  if(id==2)density=1e11_wp
  ns=density;cv=density*kb/(gammas(6)-1)
  do window=3,3
    if(window==2.and.id==2)cycle
    total=.0001_wp*(density/1e5_wp)
    if(window==2)total=.8_wp
    do method=3,3
      do level=0,2
        nt=8*2**level;dt=total/nt
        write(name,'(a,i0,a,i0,a,i0,a,i0,a)') 'conduction_d',id,'_w',window,'_m',method,'_l',level,'.csv'
        open(newunit=u,file=trim(directory)//'/'//trim(name),status='new')
        write(u,'(a)') 'step,node,t,dt,ns,dx,Cv,Told,Traw,Tnew,lambda_applied,reaction,drift,conduction,source,reservoir'
        temp=600
        do i=1,n
          temp(i,1,1)=600+900._wp*(i-1)/(n-1)
        enddo
        do it=1,nt
          old=temp
          raw=midpoint_TRBDF23D(6,x,ns,temp,nn,j1,tn,1500._wp,dt,inc,lam)
          temp=max(raw,100._wp)
          do i=1,n
            write(u,'(i0,",",i0,14(",",es25.17))') it,i,it*dt,dt,density,dx,cv,old(i,1,1), &
              raw(i,1,1),temp(i,1,1),lam(i,1,1),inc(i,1,1,:)
          enddo
        enddo
        close(u)
      enddo
    enddo
  enddo
enddo
! The manually allocated probe mesh does not claim a complete application grid.
deallocate(x%h1,x%h2,x%h3,x%dx1,x%dx1i)
end subroutine

subroutine kirchhoff_flux()
integer :: level,n,i,u
real(wp),allocatable :: temp(:,:,:),ns(:,:,:),nn(:,:,:,:),j1(:,:,:),lam(:,:,:),beta(:,:,:)
real(wp) :: dx,xx
open(newunit=u,file=trim(directory)//'/kirchhoff_flux.csv',status='new')
write(u,'(a)') 'level,node,n,x,dx,T,lambda'
do level=0,2
  n=64*2**level+1;dx=1e5_wp/(n-1)
  allocate(temp(-1:n+2,-1:3,-1:3),ns(-1:n+2,-1:3,-1:3),j1(-1:n+2,-1:3,-1:3))
  allocate(nn(n,1,1,ln),lam(n,1,1),beta(n,1,1))
  temp=600;ns=1e11_wp;j1=0;nn=0
  do i=1,n
    xx=real(i-1,wp)/real(n-1,wp)
    temp(i,1,1)=(600._wp**3.5_wp+(1500._wp**3.5_wp-600._wp**3.5_wp)*xx)**(1/3.5_wp)
  enddo
  call thermal_conduct(6,temp,ns,nn,j1,lam,beta)
  do i=1,n
    write(u,'(i0,",",i0,",",i0,4(",",es25.17))') level,i,n,(i-1)*dx,dx,temp(i,1,1),lam(i,1,1)
  enddo
  deallocate(temp,ns,j1,nn,lam,beta)
enddo
close(u)
end subroutine
subroutine wrapper_checks()
integer,parameter :: n=49
type(cartmesh) :: x
real(wp) :: ns(-1:n+2,-1:3,-1:3,lsp),temp(-1:n+2,-1:3,-1:3,lsp)
real(wp) :: old(-1:n+2,-1:3,-1:3),raw(-1:n+2,-1:3,-1:3),plain(-1:n+2,-1:3,-1:3)
real(wp) :: nn(n,1,1,ln),j(-1:n+2,-1:3,-1:3),tn(n,1,1),inc(n,1,1,5)
integer :: i,flag
call set_total_grid_sizes(n,1,1);call set_subgrid_sizes(1,1)
x%lx1=n;x%lx2=1;x%lx3=1
allocate(x%h1(-1:n+2,-1:3,-1:3),x%h2(-1:n+2,-1:3,-1:3),x%h3(-1:n+2,-1:3,-1:3))
allocate(x%dx1(0:n+2),x%dx1i(n))
x%h1=1;x%h2=1;x%h3=1;x%dx1=1e5_wp/(n-1);x%dx1i=x%dx1(1)
ns=1e5_wp;nn=0;j=0;tn=600;temp=600
call set_gridflag(2)
do i=1,n
 temp(i,1,1,:)=600+900._wp*(i-1)/(n-1)
enddo
old=temp(:,:,:,6)
raw=midpoint_TRBDF23D(6,x,ns(:,:,:,6),old,nn,j,tn,1500._wp,1.25e-5_wp,inc)
plain=midpoint_TRBDF23D(6,x,ns(:,:,:,6),old,nn,j,tn,1500._wp,1.25e-5_wp)
if(any(raw/=plain))error stop 'optional midpoint diagnostics changed the state'
if(midpoint_conduction_enabled) then
 call energy_diffusion(1.25e-5_wp,x,ns,temp,j,nn,tn,2,1500._wp)
 if(any(temp(:,:,:,6)/=max(raw,100._wp)))error stop 'standard production path differs from tested midpoint operator'
endif
! Inverted and standard zero-flux ion boundaries on constant temperature.
do flag=1,2
 call set_gridflag(flag);old=600
 raw=midpoint_TRBDF23D(6,x,ns(:,:,:,6),old,nn,j,tn,-1e-6_wp,1.25e-5_wp,inc)
 if(maxval(abs(raw(1:n,1,1)-600._wp))>1e-11_wp)error stop 'midpoint Neumann constant state changed'
enddo
deallocate(x%h1,x%h2,x%h3,x%dx1,x%dx1i)
end subroutine wrapper_checks
end program
