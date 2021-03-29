program newton_testdriver

use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use dipole, only : qp2rtheta

implicit none

real(wp) :: q,p
real(wp) :: r=0,theta=0
real(wp) :: rcorrect=0,thetacorrect=0


! test the conversion of a single point
q=-0.35615504
p=1.2790947
rcorrect=6803179.761800971
thetacorrect=1.9891332403471418
call qp2rtheta(q,p,r,theta)
print*, 'solution:  r=',r,'; theta=',theta
print*, 'known good solution:  r=',rcorrect,'; theta=',thetacorrect

if (r-rcorrect>1e4_wp .or. theta-thetacorrect>1e-4_wp) then
  error stop ' Excessive error in coordinate conversion!'
end if

end program newton_testdriver
