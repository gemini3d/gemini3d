program audit_coordinates
use phys_consts, only: wp,pi
use geomagnetic, only: geog2geomag,geomag2geog,rotgg2gm,rotgm2gg,ECEFspher2ENU
implicit none
real(wp) :: theta1,theta2,phi1,phi2,lon,lat,ident(3,3)
real(wp) :: alt(2,2,2),theta(2,2,2),phi(2,2,2),x(2,2,2),y(2,2,2),z(2,2,2)
real(wp) :: xsmall(2,1,2),ysmall(2,1,2),zsmall(2,1,2)
integer :: i,j
! Equivalent longitudes in three conventions must produce the same result.
do i=-720,720,15
 call geog2geomag(real(i,wp),65._wp,phi1,theta1)
 call geog2geomag(modulo(real(i,wp),360._wp),65._wp,phi2,theta2)
 if(abs(phi1-phi2)>1e-12_wp.or.abs(theta1-theta2)>1e-12_wp)error stop 'longitude equivalence'
 call geomag2geog(phi1,theta1,lon,lat)
 if(abs(lat-65._wp)>1e-10_wp.or.abs(modulo(lon-real(i,wp)+180,360._wp)-180)>1e-10_wp) &
   error stop 'coordinate round trip'
enddo
ident=matmul(rotgg2gm(),rotgm2gg())
do i=1,3
 do j=1,3
  if(abs(ident(i,j)-merge(1._wp,0._wp,i==j))>1e-14_wp)error stop 'rotation inverse'
 enddo
enddo
! A 3D call must not leave SAVE state that changes subsequent collapsed-dimension calls.
alt=100000;theta=pi/3;phi=0.2_wp
call ECEFspher2ENU(alt,theta,phi,pi/3,0._wp,x,y,z)
call ECEFspher2ENU(alt(:,1:1,:),theta(:,1:1,:),phi(:,1:1,:),pi/3,0._wp,xsmall,ysmall,zsmall)
if(maxval(abs(xsmall))>1e-8_wp)error stop '3D state leaked into 2D transform'
print *,'coordinate checks passed'
end program
