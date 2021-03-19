program newton_testdriver

use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use dipole, only : qp2rtheta

implicit none

integer, parameter :: lq=512,lp=512
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(2) :: qlims=[-0.5851937,0.5851937]
real(wp), dimension(2) :: plims=[1.2053761,1.5820779]
real(wp), dimension(lq,lp) :: r,theta
!real(wp) :: rcorrect=0,thetacorrect=0
integer :: iq,ip


! define a grid
q=[(qlims(1) + (qlims(2)-qlims(1)/lq*(iq-1)),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1)/lp*(ip-1)),ip=1,lp)]

!print*, q
!print*, p


! test the conversion of a single point
!q=-0.35615504
!p=1.2790947
!rcorrect=6803179.761800971
!thetacorrect=1.9891332403471418
do iq=1,lq
  do ip=1,lp
    call qp2rtheta(q(ip),p(ip),r(iq,ip),theta(iq,ip))
    !print*, 'solution:  r=',r(iq,ip),'; theta=',theta(iq,ip)
    !print*, 'matlab solution:  r=',rcorrect,'; theta=',thetacorrect
  end do
end do

end program newton_testdriver
