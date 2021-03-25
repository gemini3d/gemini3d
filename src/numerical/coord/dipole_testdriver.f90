program newton_testdriver

use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use dipole, only : qp2rtheta,rtheta2qp

implicit none

integer, parameter :: lq=16,lp=16
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(2) :: qlims=[-0.5851937,0.5851937]
real(wp), dimension(2) :: plims=[1.2053761,1.5820779]
real(wp), dimension(lq,lp) :: r,theta,qtest,ptest,errq,errp
real(wp) :: maxerrq,maxerrp
integer :: iq,ip


! define a grid
q=[(qlims(1) + (qlims(2)-qlims(1)/lq*(iq-1)),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1)/lp*(ip-1)),ip=1,lp)]


! test the conversion of a set of grid points
do iq=1,lq
  do ip=1,lp
    call qp2rtheta(q(iq),p(ip),r(iq,ip),theta(iq,ip))
  end do
end do

! now convert back to dipole coords so we can check
do iq=1,lq
  do ip=1,lp
    call rtheta2qp(r(iq,ip),theta(iq,ip),qtest(iq,ip),ptest(iq,ip))
  end do
end do

! now compute errors for each conversion that has been done
do iq=1,lq
  do ip=1,lp
    errq(iq,ip)=qtest(iq,ip)-q(iq)
    errp(iq,ip)=ptest(iq,ip)-p(ip)
  end do
end do

! provide some indication to the call program as to whether this test passed or not
maxerrq=maxval(abs(errq))
maxerrp=maxval(abs(errp))
print*, ' Maximum error in computed q and p:  ',maxerrq,maxerrp
if (maxerrq>1e-4_wp .or. maxerrp>1e-4_wp) error stop ' Excessive error in coordinate conversion!'

end program newton_testdriver
