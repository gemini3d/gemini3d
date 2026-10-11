program diffusion_scale
use phys_consts, only: wp
use PDEparabolic, only: backEuler1D,TRBDF21D
implicit none(type,external)
integer,parameter::n=49
real(wp)::u(n),initial(n),a(n),b(n),c(n),d(n),e(n),dx(0:n+2),dxi(n),dt,err,bcerr
integer::i,method
logical::passed
a=0;b=0;c=1e18_wp;d=3e-5_wp;e=0
dx=1e5_wp/(n-1);dxi=dx(1);dt=6.25e-5_wp
initial=[(600+900._wp*(i-1)/(n-1),i=1,n)]
passed=.true.
do method=1,2
 u=initial
 if(method==1)then
   u=backEuler1D(u,a,b,c,d,e,600._wp,1500._wp,dt,[0,0],dx,dxi)
 else
   u=TRBDF21D(u,a,b,c,d,e,600._wp,1500._wp,dt,[0,0],dx,dxi)
 endif
 err=maxval(abs(u-initial));bcerr=max(abs(u(1)-600),abs(u(n)-1500))
 print '(i0,2(1x,es25.17))',method,err,bcerr
 passed=passed.and.bcerr<=1e-11_wp
enddo
if(.not.passed)error stop 'prescribed boundary temperature lost under SI coefficient scaling'
end program
