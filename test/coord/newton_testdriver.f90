program newton_testdriver

use phys_consts, only: wp
use meshobj_dipole, only : qp2rtheta

implicit none (type, external)

real(wp),dimension(1) :: q,p,phi
real(wp) :: r=0,theta=0
real(wp) :: rcorrect=0,thetacorrect=0

! test the conversion of a single point
q(1)=-0.35615504
p(1)=1.2790947
phi(1)=3.14
rcorrect=6803179.761800971_wp
thetacorrect=1.9891332403471418_wp
call qp2rtheta(q(1),p(1),r,theta)
print*, 'solution:  r=',r,'; theta=',theta
print*, 'known good solution:  r=',rcorrect,'; theta=',thetacorrect

if (r-rcorrect>1e4_wp .or. theta-thetacorrect>1e-4_wp) then
  error stop ' Excessive error in coordinate conversion!'
end if

end program newton_testdriver
