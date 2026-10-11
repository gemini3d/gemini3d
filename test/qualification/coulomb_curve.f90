program coulomb_curve
use phys_consts,only:wp,lsp,kb,ms
use collisions,only:coulomb_colln
implicit none(type,external)
real(wp),parameter::values(13)=[0._wp,1e-8_wp,1e-6_wp,.001_wp,.01_wp,.05_wp, &
 .09999999_wp,.1_wp,.10000001_wp,.5_wp,1._wp,3._wp,10._wp]
real(wp)::ns(-1:3,-1:3,-1:3,lsp),ts(-1:3,-1:3,-1:3,lsp),v(-1:3,-1:3,-1:3,lsp)
real(wp)::nu(1,1,1),phi(1,1,1),psi(1,1,1),speed,w,oracle,error
integer::i
ns=0;ts=1000;v=0;ns(:,:,:,1)=2e11_wp;ns(:,:,:,6)=8e10_wp
speed=sqrt(2*kb*1000/(ms(1)*ms(6)/(ms(1)+ms(6))))
do i=1,size(values)
 v(:,:,:,1)=values(i)*speed;w=abs(v(1,1,1,1))/speed
 call coulomb_colln(1,6,ns,ts,v,nu,phi,psi)
 ! Independent numerical integral on a fixed interval. This does not reuse
 ! either production power-series coefficients or the erf expression.
 oracle=integral(w);error=abs(phi(1,1,1)-oracle)
 print '(2(es25.17,1x))',w,error
 if(error>5e-14_wp)error stop 'Coulomb drag integral and continuity'
enddo
contains
pure function integral(w)result(value)
real(wp),intent(in)::w
real(wp)::value,total,compensation,term,y,next,x,weight
integer,parameter::n=4096
integer::j
total=exp(-w*w);compensation=0
do j=1,n-1
 x=real(j,wp)/n;weight=2
 if(modulo(j,2)==1)weight=4
 term=weight*x*x*exp(-(w*x)**2)
 y=term-compensation;next=total+y
 compensation=(next-total)-y;total=next
enddo
value=total/n
end function
end program
