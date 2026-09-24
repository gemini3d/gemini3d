program transport_convergence
use phys_consts, only: wp,pi
use advec, only: advec1D_MC_curv
implicit none(type,external)
integer :: n,level,i,it,nt
real(wp) :: h,dt,error(3),order,flux_residual,max_residual
real(wp), allocatable :: f(:),old(:),v(:),dx(:),dxi(:),hc(:),hi(:),flux(:),exact(:)
do level=1,3
  n=32*2**(level-1);nt=2*n;h=1._wp/n;dt=1._wp/nt
  allocate(f(-1:n+2),old(-1:n+2),v(n+1),dx(0:n+2),dxi(n),hc(-1:n+2),hi(n+1),flux(n+1),exact(n))
  v=1;dx=h;dxi=h;hc=1;hi=1
  do i=1,n
    f(i)=2+sin(2*pi*(i-.5_wp)*h)
  enddo
  exact=f(1:n);max_residual=0
  do it=1,nt
    f(-1:0)=f(n-1:n);f(n+1:n+2)=f(1:2);old=f
    f=advec1D_MC_curv(old,v,dt,dx,dxi,hc,hi,hi,flux)
    flux_residual=abs(sum(f(1:n)-old(1:n))*h+flux(n+1)-flux(1))/sum(abs(old(1:n))*h)
    max_residual=max(max_residual,flux_residual)
  enddo
  error(level)=sum(abs(f(1:n)-exact))*h
  print '(A,I0,2(A,ES24.16))','cells=',n,' L1=',error(level),' balance=',max_residual
  if(max_residual>1e-11_wp) error stop 'Periodic transport balance budget exceeded'
  deallocate(f,old,v,dx,dxi,hc,hi,flux,exact)
enddo
do level=2,3
  order=log(error(level-1)/error(level))/log(2._wp)
  print '(A,F10.6)','measured advection space/time order: ',order
  if(order<1.8_wp) error stop 'Smooth advection convergence below frozen budget'
enddo
! Curved weights and known nonzero boundary flux (no inference from state delta).
n=12
allocate(f(-1:n+2),old(-1:n+2),v(n+1),dx(0:n+2),dxi(n),hc(-1:n+2),hi(n+1),flux(n+1))
f=2;v=1;dx=.1_wp;dxi=.1_wp;hc=2
hi=[(1._wp+.01_wp*i,i=1,n+1)]
old=f;dt=.001_wp
f=advec1D_MC_curv(old,v,dt,dx,dxi,hc,hi,hi,flux)
flux_residual=abs(sum((f(1:n)-old(1:n))*hc(1:n)*dxi)+2*dt*(hi(n+1)-hi(1)))
if(flux_residual>1e-11_wp) error stop 'Known curved-boundary flux budget exceeded'
print '(A,ES24.16)','curved known-boundary residual: ',flux_residual
end program
