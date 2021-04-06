program grid_testdriver

use, intrinsic :: ISO_Fortran_env,  only : wp=>real64
use meshobj_dipole, only : dipolemesh

implicit none

integer, parameter :: lq=1024,lp=1024
real(wp), dimension(lq) :: q
real(wp), dimension(lp) :: p
real(wp), dimension(2) :: qlims=[-0.5851937,0.5851937]
real(wp), dimension(2) :: plims=[1.2053761,1.5820779]
real(wp), dimension(1) :: phi=[3.14]

real(wp), dimension(lq,lp) :: r,theta,qtest,ptest,errq,errp
real(wp) :: maxerrq,maxerrp
integer :: iq,ip,irepeat
type(dipolemesh) :: x     ! dummy object to allow us to access methods for coodinate conversion


! define a grid
q=[(qlims(1) + (qlims(2)-qlims(1))/lq*(iq-1),iq=1,lq)]
p=[(plims(1) + (plims(2)-plims(1))/lp*(ip-1),ip=1,lp)]

! initialize object so we can call its type-bound procedures
call x%set_coords(q,p,phi,p,phi)
call x%init_dipolemesh()

! test the conversion of a set of grid points
do irepeat=1,3    ! do mulitple times as would be require for a full grid generation (to compare times against scripts)
  call x%calc_rtheta_2D(q,p,r,theta)
end do

! now convert back to dipole coords so we can check
call x%calc_qp_2D(r,theta,qtest,ptest)

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
print*, ' Test system size:  ',lq,lp
print*, ' Maximum error in computed q and p:  ',maxerrq,maxerrp
if (maxerrq>1e-6_wp .or. maxerrp>1e-6_wp) error stop ' Excessive error in coordinate conversion!'

end program grid_testdriver
