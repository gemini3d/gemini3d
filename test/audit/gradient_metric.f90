program qualification_gradient_metric
use phys_consts, only: wp
use meshobj_cart, only: cartmesh
use calculus, only: grad3D2
implicit none
type(cartmesh) :: x
real(wp) :: f(2,4,3),g(2,4,3),expected(2,4,3)
integer :: j,k
x%lx1=2;x%lx2=4;x%lx3=3;x%lx2all=4;x%lx3all=3
allocate(x%h2(2,4,3),x%dx2(4));x%dx2=2._wp
! Manufactured scalar linear in coordinate x2, with a metric that varies
! across x3. Its physical derivative is 1/h2 at every boundary and interior.
do k=1,3
  x%h2(:,:,k)=real(k,wp)
  do j=1,4
    f(:,j,k)=2._wp*real(j,wp)
  enddo
enddo
expected=1._wp/x%h2
g=grad3D2(f,x,1,2,1,4,1,3)
if(maxval(abs(g-expected))>1e-14_wp) error stop 'Incorrect metric at x2 boundary'
deallocate(x%h2,x%dx2)
print *, 'Manufactured curved-metric gradient passed'
end program
