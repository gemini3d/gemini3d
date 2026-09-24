! Audit addition 2026-09-16. Apache-2.0.
program audit_interp
use phys_consts, only: wp
use interpolation, only: interp1,interp2,interp3
implicit none
real(wp) :: x(2),q(4),a(2),b(2,2),c(2,2,2),y(4),expected(4)
integer :: i,j,k
x=[0._wp,1._wp];q=[0._wp,0.5_wp,1._wp,2._wp]
a=2*x+3
expected=[3._wp,4._wp,5._wp,0._wp]
y=interp1(x,a,q)
if (maxval(abs(y-expected))>1e-12_wp) error stop 'interp1 affine/boundary/OOD'
do j=1,2
 do i=1,2
  b(i,j)=2*x(i)+3*x(j)+1
 enddo
enddo
y=interp2(x,x,b,q,q);expected=[1._wp,3.5_wp,6._wp,0._wp]
if (maxval(abs(y-expected))>1e-12_wp) error stop 'interp2 affine/boundary/OOD'
do k=1,2
 do j=1,2
  do i=1,2
   c(i,j,k)=2*x(i)+3*x(j)+4*x(k)+1
  enddo
 enddo
enddo
y=interp3(x,x,x,c,q,q,q);expected=[1._wp,5.5_wp,10._wp,0._wp]
if (maxval(abs(y-expected))>1e-12_wp) error stop 'interp3 affine/boundary/OOD'
! A singleton 1D axis represents an invariant driver dimension.
y=interp1([0._wp],[7._wp],q)
if (any(y/=7._wp)) error stop 'interp1 singleton invariant dimension'
print *, 'PASS affine interpolation on 2-point axes, endpoints and zero-fill outside domain'
end program
