program energy_operator_balance
use phys_consts, only: wp
use PDEparabolic, only: backEuler1D,TRBDF21D
use diffusion, only: backEuler3D,TRBDF23D
use meshobj_cart, only: cartmesh
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
integer,parameter :: n=13
real(wp) :: old(n),y(n),a(n),b(n),c(n),d(n),e(n),dx(0:n+2),dxi(n),terms(n,5),dt,err,scale
real(wp) :: old3(-1:n+2,-1:3,-1:3),new3(-1:n+2,-1:3,-1:3)
real(wp) :: aa(n,1,1),zz(n,1,1),cc(n,1,1),terms3(n,1,1,5)
type(cartmesh) :: x
integer :: i,method,bc1,bc2,probe
! Nonuniform mesh, variable coefficients, and both boundary types. Independently
! reconstruct the discrete PDE budget, including prescribed nodal endpoints.
dx=[(0.7_wp+0.03_wp*i,i=0,n+2)];dxi=0.5_wp*(dx(1:n)+dx(2:n+1))
do probe=1,3
 if(probe==1)dt=0.0001_wp
 if(probe==2)dt=0.1_wp
 if(probe==3)dt=2._wp
 old=[(2+sin(real(i,wp)),i=1,n)]
 a=-0.3_wp;b=0.12_wp;c=0.8_wp;d=[(1+0.02_wp*i,i=1,n)];e=0.9_wp
 do method=1,2
  do bc1=0,1
   do bc2=0,1
    if(method==1)then
      y=backEuler1D(old,a,b,c,d,e,1.2_wp,2.4_wp,dt,[bc1,bc2],dx,dxi,increments=terms)
    else
      y=TRBDF21D(old,a,b,c,d,e,1.2_wp,2.4_wp,dt,[bc1,bc2],dx,dxi, &
        Tsmin_mid=0.6_wp,Tsmax_mid=1.7_wp,increments=terms)
    endif
    scale=max(1._wp,maxval(abs(y)),maxval(abs(old)))
    err=maxval(abs(y-old-sum(terms,2)))/scale
    if(.not.all(ieee_is_finite(y)).or.err>1e-11_wp) error stop 'Parabolic stage balance'
    if(any(terms(2:n-1,5)/=0)) error stop 'Interior residual attributed to boundary'
    if(any(terms(1,1:4)/=0).or.any(terms(n,1:4)/=0)) error stop 'Boundary solved as PDE'
   enddo
  enddo
 enddo
enddo
! Source-free pure conduction must telescope with reciprocal C volume weights.
a=0;b=0;e=0;dt=0.1_wp
c=[(0.5_wp+0.03_wp*i,i=1,n)]
y=backEuler1D(old,a,b,c,d,e,1.2_wp,2.4_wp,dt,[0,0],dx,dxi,increments=terms)
err=abs(sum(terms(2:n-1,3)*dxi(2:n-1)/c(2:n-1)) - dt* &
 (0.5_wp*(d(n)+d(n-1))*(y(n)-y(n-1))/dx(n) - 0.5_wp*(d(1)+d(2))*(y(2)-y(1))/dx(2)))
if(err>1e-11_wp) error stop 'Conductive face telescoping'
! Wrappers must preserve every halo cell; compare requested and absent ledger.
allocate(x%dx1(0:n+2),x%dx1i(n))
x%dx1=dx;x%dx1i=dxi
old3=321._wp;old3(1:n,1,1)=old;old3(0,1,1)=1.2_wp;old3(n+1,1,1)=2.4_wp
aa=-0.1_wp;zz=0;cc=1
new3=TRBDF23D(old3,aa,zz,cc,cc,zz,dt,x,terms3)
if(any(new3(-1,:,:)/=old3(-1,:,:)).or.any(new3(:,0,:)/=old3(:,0,:))) error stop 'TR halos undefined'
err=maxval(abs(new3(1:n,1,1)-old-sum(terms3(:,1,1,:),2)))
if(err>1e-11_wp) error stop '3D TR ledger'
new3=backEuler3D(old3,aa,zz,cc,cc,zz,dt,x)
if(any(new3(n+2,:,:)/=old3(n+2,:,:)).or.any(new3(:,:,3)/=old3(:,:,3))) error stop 'BE halos undefined'
print '(a)', 'PASS: 24 parabolic cases; face telescoping; 3D ledgers and defined halos'
end program
