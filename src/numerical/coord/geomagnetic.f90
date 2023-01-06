module geomagnetic

!> transformations and data relevant to converting geographic to geomagnetic coordinates

use phys_consts, only: wp,Re,pi
implicit none (type, external)


! magnetic pole location in geographic coordinates
real(wp), parameter :: thetan=11*pi/180
real(wp), parameter :: phin=289*pi/180

private
public :: geomag2geog, geog2geomag, r2alt, alt2r, rotgg2gm, rotgm2gg, ECEFspher2ENU

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


  !> take a set of ECEF spherical coordinatese and convert into ENU (geog or geom)
  subroutine ECEFspher2ENU(alt,theta,phi,theta1,phi1,x,y,z)
    real(wp), intent(in), dimension(:,:,:) :: alt,theta,phi
    real(wp), intent(in) :: theta1,phi1
    real(wp), intent(inout), dimension(:,:,:) :: x,y,z    !ENU
    real(wp) :: xp,yp
    real(wp) :: theta2,theta3,gamma1,gamma2,phi2,phi3
    integer :: lx1,lx2,lx3,ix1,ix2,ix3    ! local copies

    lx1=size(alt,1)
    lx2=size(alt,2)
    lx3=size(alt,3)
    if (size(z,1)/=lx1 .or. size(z,2)/=lx2 .or. size(z,3)/=lx3) error stop 'ECEFspher2ENU:  inconsistent input array sizes'

    z(:,:,:)=alt(:,:,:)
    do ix3=1,lx3
      do ix2=1,lx2
        do ix1=1,lx1
          theta2=theta(ix1,ix2,ix3)                    !field point zenith angle
      
          if (size(alt,2)/=1) then
            phi2=phi(ix1,ix2,ix3)                      !field point azimuth, full 3D calculation
          else
            phi2=phi1                                    !assume the longitude is the samem as the source in 2D, i.e. assume the source epicenter is in the meridian of the grid
          end if
      
          !we need a phi locationi (not spherical phi, but azimuth angle from epicenter), as well, but not for interpolation - just for doing vector rotations
          theta3=theta2
          phi3=phi1
          gamma1=cos(theta2)*cos(theta3)+sin(theta2)*sin(theta3)*cos(phi2-phi3)
          if (gamma1 > 1) then     !handles weird precision issues in 2D
            gamma1 = 1
          else if (gamma1 < -1) then
            gamma1 = -1
          end if
          gamma1=acos(gamma1)
      
          gamma2=cos(theta1)*cos(theta3)+sin(theta1)*sin(theta3)*cos(phi1-phi3)
          if (gamma2 > 1) then     !handles weird precision issues in 2D
            gamma2= 1
          else if (gamma2 < -1) then
            gamma2= -1
          end if
          gamma2=acos(gamma2)
          xp=Re*gamma1
          yp=Re*gamma2     !this will likely always be positive, since we are using center of earth as our origin, so this should be interpreted as distance as opposed to displacement
      
          ! coordinates from distances
          if (theta3>theta1) then       !place distances in correct quadrant, here field point (theta3=theta2) is is SOUTHward of source point (theta1), whreas yp is distance northward so throw in a negative sign
            yp=-yp            !do we want an abs here to be safe
          end if
          if (phi2<phi3) then     !assume we aren't doing a global grid otherwise need to check for wrapping, here field point (phi2) less than source point (phi3=phi1)
            xp=-xp
          end if

          x(ix1,ix2,ix3)=xp
          y(ix1,ix2,ix3)=yp
        end do
      end do
    end do
  end subroutine ECEFspher2ENU
end module geomagnetic
