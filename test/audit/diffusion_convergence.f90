program qualification_diffusion
use phys_consts, only: wp,pi
use PDEparabolic, only: backEuler1D,TRBDF21D
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
real(wp) :: spatial(3),temporal(3,2),order,err
integer :: level,method
! Thresholds are declared from method order before running: spatial >=1.8,
! backward Euler time >=0.9, TRBDF2 time >=1.8. These qualify this operator.
do level=1,3
  call solve(2**(level+3)+1,1000,0.001_wp,2,.true.,spatial(level))
  do method=1,2
    call solve(9,2**(level+3),1._wp,method,.false.,temporal(level,method))
  enddo
enddo
print '(a,3es24.16)', 'SPACE_ERRORS ',spatial
do level=1,2
  order=log(spatial(level)/spatial(level+1))/log(2._wp)
  if (.not.ieee_is_finite(order).or.order<1.8_wp) error stop 'diffusion spatial order'
enddo
do method=1,2
  print '(a,i0,3es24.16)', 'TIME_ERRORS ',method,temporal(:,method)
  do level=1,2
    order=log(temporal(level,method)/temporal(level+1,method))/log(2._wp)
    if (.not.ieee_is_finite(order)) error stop 'nonfinite temporal convergence'
    if (method==1.and.order<0.9_wp) error stop 'backward Euler temporal order'
    if (method==2.and.order<1.8_wp) error stop 'TRBDF2 temporal order'
  enddo
enddo
! Manufactured constant with a volumetric source and matching boundary values.
call source_balance(err)
if(err>1e-11_wp)error stop 'source and boundary balance'
print '(a,es24.16)', 'SOURCE_BALANCE_ERROR ',err
contains
subroutine solve(n,nt,total,method,spatial_case,error)
integer,intent(in)::n,nt,method
real(wp),intent(in)::total
logical,intent(in)::spatial_case
real(wp),intent(out)::error
real(wp)::u(n),xx(n),a(n),b(n),c(n),d(n),e(n),dx(0:n+2),dxi(n),exact(n),dt,bc
integer::i,it
xx=[(real(i-1,wp)/real(n-1,wp),i=1,n)]
dx=1._wp/(n-1);dxi=dx(1);dt=total/nt
b=0;c=1;e=0
if(spatial_case)then
  u=sin(2*pi*xx);a=0;d=1;exact=exp(-4*pi*pi*total)*u
else
  u=1;a=-1;d=0;exact=exp(-total)
endif
do it=1,nt
  bc=0
  if(.not.spatial_case)bc=exp(-it*dt)
  if(method==1)then
    u=backEuler1D(u,a,b,c,d,e,bc,bc,dt,[0,0],dx,dxi)
  else
    u=TRBDF21D(u,a,b,c,d,e,bc,bc,dt,[0,0],dx,dxi)
  endif
enddo
error=sqrt(sum((u-exact)**2)/sum(exact**2))
if(.not.ieee_is_finite(error).or.error<=0)error stop 'invalid refinement error'
end subroutine
subroutine source_balance(error)
real(wp),intent(out)::error
real(wp)::u(9),z(9),one(9),dx(0:11),dxi(9),e(9),dt
integer::i
u=2;z=0;one=1;dx=0.125_wp;dxi=dx(1);e=3;dt=0.01_wp
do i=1,100
  u=TRBDF21D(u,z,z,one,one,e,2+3*i*dt,2+3*i*dt,dt,[0,0],dx,dxi, &
    Tsmin_mid=2+3*(i-0.5_wp)*dt,Tsmax_mid=2+3*(i-0.5_wp)*dt)
enddo
error=maxval(abs(u-5._wp))
end subroutine
end program
