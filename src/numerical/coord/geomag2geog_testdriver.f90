program geomag2geog_testdriver

use phys_consts, only: wp,Re,pi
use geomagnetic, only : geomag2geog,geog2geomag,r2alt,alt2r

implicit none (type, external)

integer, parameter :: lr=128,ltheta=192,lphi=96
integer :: ir,itheta,iphi
real(wp), dimension(2) :: rlims, thetalims, philims
real(wp), dimension(lr) :: rarr
real(wp), dimension(ltheta) :: thetaarr
real(wp), dimension(lphi) :: phiarr
real(wp), dimension(:,:,:), allocatable :: r,theta,phi,  glon,glat,alt, rtest,thetatest,phitest
real(wp) :: r0,theta0,phi0,glon0,glat0,alt0
real(wp) :: r0test,theta0test,phi0test
real(wp) :: maxerrr,maxerrtheta,maxerrphi

allocate(r(lr,ltheta,lphi))
allocate(theta, phi, glon, glat, alt, rtest, thetatest, phitest, mold=r)

!> define points of interest
r0=Re+80e3
theta0=pi/36
phi0=pi/36
rlims=[r0,r0+900e3]
thetalims=[theta0,pi-theta0]
philims=[phi0,2*pi-phi0]
rarr=[(r0 + (rlims(2)-rlims(1))/lr*(ir-1),ir=1,lr)]
thetaarr=[(theta0 + (thetalims(2)-thetalims(1))/ltheta*(itheta-1),itheta=1,ltheta)]
phiarr=[(phi0 + (philims(2)-philims(1))/lphi*(iphi-1),iphi=1,lphi)]
do ir=1,lr
  do itheta=1,ltheta
    do iphi=1,lphi
      r(ir,itheta,iphi)=rarr(ir)
      theta(ir,itheta,iphi)=thetaarr(itheta)
      phi(ir,itheta,iphi)=phiarr(iphi)
    end do
  end do
end do

!> test the conversion of a single point
call geomag2geog(phi0,theta0,glon0,glat0)
alt0=r2alt(r0)
call geog2geomag(glon0,glat0,phi0test,theta0test)
r0test=alt2r(alt0)
print*, ' r,theta,phi coordinate diff:  ',abs(r0-r0test),abs(theta0-theta0test),abs(phi0-phi0test)
if (abs(phi0-phi0test)>1e-6 .or. abs(theta0-theta0test)>1e-6 .or. abs(r0-r0test)>1e-6) then
  error stop ' excessive error in single point coordinate conversion'
end if

!> test a 3D array conversion
call geomag2geog(phi,theta,glon,glat)
alt=r2alt(r)
print*, minval(glon),maxval(glon),minval(glat),maxval(glat)
call geog2geomag(glon,glat,phitest,thetatest)
rtest=alt2r(alt)
maxerrr=maxval(abs(rtest-r))
maxerrtheta=maxval(abs(thetatest-theta))
maxerrphi=maxval(abs(phitest-phi))
print*, ' r,theta,phi array max diff:  ',maxerrr,maxerrtheta,maxerrphi
if ( any([maxerrr,maxerrtheta,maxerrphi]>1e-6) ) error stop ' Excessive error in array coordinate conversion!'

end program geomag2geog_testdriver
