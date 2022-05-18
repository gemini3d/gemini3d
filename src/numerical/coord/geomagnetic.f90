module geomagnetic

!> transformations and data relevant to converting geographic to geomagnetic coordinates

use phys_consts, only: wp,Re,pi
implicit none (type, external)


! magnetic pole location in geographic coordinates
real(wp), parameter :: thetan=11*pi/180
real(wp), parameter :: phin=289*pi/180

private
public :: geomag2geog, geog2geomag, r2alt, alt2r, rotgg2gm, rotgm2gg

contains
  !> convert geomagnetic coordinates to geographic
  elemental subroutine geog2geomag(glon,glat,phi,theta)
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
  end subroutine geog2geomag
  
  
  !> convert geographic coordinates to geomagnetic; do not use at the magnetic pole!!!
  elemental subroutine geomag2geog(phi,theta,glon,glat)
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
  end subroutine geomag2geog
  
  
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
  
  
  !> return elemental rotation matrix for rotations about z-axis
  function rotz(alpha) result(Rz)
    real(wp), intent(in) :: alpha     ! must be in radians
    real(wp), dimension(3,3) :: Rz
  
    Rz(1:3,1:3)=0._wp
    Rz(1,1)=cos(alpha)
    Rz(1,2)=-sin(alpha)
    Rz(2,1)=sin(alpha)
    Rz(2,2)=cos(alpha)
    Rz(3,3)=1._wp
  end function rotz
  
  
  !> return elemental rotation matrix for rotations about y-axis
  function roty(alpha) result(Ry)
    real(wp), intent(in) :: alpha     ! must be in radians
    real(wp), dimension(3,3) :: Ry
  
    Ry(1:3,1:3)=0._wp
    Ry(1,1)=cos(alpha)
    Ry(1,3)=sin(alpha)
    Ry(2,2)=1._wp
    Ry(3,1)=-sin(alpha)
    Ry(3,3)=cos(alpha)
  end function roty
  
  
  !> a rotation matrix to go from geomagnetic ECEF to geographic ECEF
  function rotgm2gg() result(Rgm2gg)
    real(wp), dimension(3,3) :: Rgm2gg
  
    Rgm2gg=matmul(rotz(phin),roty(thetan))
  end function rotgm2gg


  !> rotation matrix to go from geographic ECEF to geomagnetic ECEF
  function rotgg2gm() result(Rgg2gm)
    real(wp), dimension(3,3) :: Rgg2gm

    Rgg2gm=matmul(transpose(roty(thetan)),transpose(rotz(phin)))
  end function rotgg2gm
end module geomagnetic
