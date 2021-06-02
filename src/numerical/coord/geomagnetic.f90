module geomagnetic

!> transformations and data relevant to converting geographic to geomagnetic coordinates

use phys_consts, only: wp,Re,pi
implicit none (type, external)


! magnetic pole location in geographic coordinates
real(wp), parameter :: thetan=11*pi/180
real(wp), parameter :: phin=289*pi/180

! overload for scalar vs. rank3 arrays
interface geomag2geog
  module procedure geomag2geog_scalar
  module procedure geomag2geog_rank3
end interface geomag2geog
interface geog2geomag
  module procedure geog2geomag_scalar
  module procedure geog2geomag_rank3
end interface geog2geomag

contains

!> convert geomagnetic coordinates to geographic
subroutine geog2geomag_rank3(glon,glat,phi,theta)
  real(wp), dimension(:,:,:), intent(in) :: glon,glat
  real(wp), dimension(:,:,:), intent(out) :: phi,theta
  real(wp), dimension(1:size(glon,1),1:size(glon,2),1:size(glon,3)) :: glonwrap
  real(wp), dimension(1:size(glon,1),1:size(glon,2),1:size(glon,3)) :: thetag
  real(wp), dimension(1:size(glon,1),1:size(glon,2),1:size(glon,3)) :: phig
  real(wp), dimension(1:size(glon,1),1:size(glon,2),1:size(glon,3)) :: argtmp,alpha

  glonwrap=mod(glon,360._wp)
  thetag=pi/2._wp-glat*pi/180._wp
  phig=glonwrap*pi/180._wp

  theta = acos(cos(thetag)*cos(thetan)+sin(thetag)*sin(thetan)*cos(phig-phin))
  argtmp = (cos(thetag)-cos(theta)*cos(thetan))/(sin(theta)*sin(thetan))
  alpha = acos( max(min(argtmp,1._wp),-1._wp) )

  where (phin>phig .and. phin-phig>pi .or. phin<phig .and. phig-phin<pi)
    phi=pi-alpha
  elsewhere
    phi=alpha+pi
  end where
end subroutine geog2geomag_rank3


subroutine geog2geomag_scalar(glon,glat,phi,theta)
  real(wp), intent(in) :: glon,glat
  real(wp), intent(out) :: phi,theta
  real(wp) :: glonwrap
  real(wp) :: thetag
  real(wp) :: phig
  real(wp) :: argtmp,alpha

  glonwrap=mod(glon,360._wp)
  thetag = pi/2 - glat*pi/180
  phig = glonwrap*pi/180

  theta = acos(cos(thetag)*cos(thetan)+sin(thetag)*sin(thetan)*cos(phig-phin))
  argtmp = (cos(thetag)-cos(theta)*cos(thetan))/(sin(theta)*sin(thetan))
  alpha = acos( max(min(argtmp,1._wp),-1._wp) )

  if (phin>phig .and. phin-phig>pi .or. phin<phig .and. phig-phin<pi) then
    phi=pi-alpha
  else
    phi=alpha+pi
  end if
end subroutine geog2geomag_scalar


!> convert geographic coordinates to geomagnetic; do not use at the magnetic pole!!!
subroutine geomag2geog_rank3(phi,theta,glon,glat)
  real(wp), dimension(:,:,:), intent(in) :: phi,theta
  real(wp), dimension(:,:,:), intent(out) :: glon,glat
  real(wp), dimension(1:size(phi,1),1:size(phi,2),1:size(phi,3)) :: thetag2p,thetag2
  real(wp), dimension(1:size(phi,1),1:size(phi,2),1:size(phi,3)) :: beta
  real(wp), dimension(1:size(phi,1),1:size(phi,2),1:size(phi,3)) :: phig2,phiwrap,argtmp

  phiwrap=mod(phi,2*pi)
  thetag2p=acos(cos(theta)*cos(thetan)-sin(theta)*sin(thetan)*cos(phiwrap))
  argtmp=(cos(theta)-cos(thetag2p)*cos(thetan))/(sin(thetag2p)*sin(thetan))
  beta=acos( max(min(argtmp,1._wp),-1._wp) )     ! deal with slight overshoots depending on precision used...

  where (phiwrap>pi)
    phig2=phin-beta
  elsewhere
    phig2=phin+beta
  end where
  phig2=mod(phig2,2*pi)
  thetag2=pi/2-thetag2p

  glon=phig2*180._wp/pi
  glat=thetag2*180._wp/pi
end subroutine geomag2geog_rank3
subroutine geomag2geog_scalar(phi,theta,glon,glat)
  real(wp), intent(in) :: phi,theta
  real(wp), intent(out) :: glon,glat
  real(wp) :: thetag2p,thetag2
  real(wp) :: beta
  real(wp) :: phig2,phiwrap
  real(wp) :: argtmp

  phiwrap=mod(phi,2*pi)
  thetag2p=acos(cos(theta)*cos(thetan)-sin(theta)*sin(thetan)*cos(phiwrap))
  argtmp=(cos(theta)-cos(thetag2p)*cos(thetan))/(sin(thetag2p)*sin(thetan))
  beta=acos( max(min(argtmp,1._wp),-1._wp) )

  if (phiwrap>pi) then
    phig2=phin-beta
  else
    phig2=phin+beta
  end if
  phig2=mod(phig2,2*pi)
  thetag2=pi/2-thetag2p

  glon=phig2*180._wp/pi
  glat=thetag2*180._wp/pi
end subroutine geomag2geog_scalar


!> convert geocentric distance into altitude (assume spherical Earth but could use other model)
elemental function r2alt(r) result(alt)
  real(wp), intent(in) :: r
  real(wp) :: alt

  alt=r-Re
end function r2alt


!> convert altitude to geocentric distance
elemental function alt2r(alt) result(r)
  real(wp), intent(in) :: alt
  real(wp) :: r

  r=alt+Re
end function alt2r

end module geomagnetic
