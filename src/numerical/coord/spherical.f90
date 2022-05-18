module spherical

!> utility routines to compute various quantities relevant to spherical coordiante systems
use phys_consts, only: wp
implicit none (type, external)

private
public :: er_spherical,etheta_spherical,ephi_spherical

!> overload calls for scalar vs. array (rank3)
interface er_spherical
  procedure er_spherical_scalar, er_spherical_rank3
end interface er_spherical
interface etheta_spherical
  procedure etheta_spherical_scalar, etheta_spherical_rank3
end interface etheta_spherical
interface ephi_spherical
  procedure ephi_spherical_scalar, ephi_spherical_rank3
end interface ephi_spherical

contains
  !> scalar:  radial unit vector (expressed in ECEF cartesian coodinates, components permuted as ix,iy,iz)
  pure function er_spherical_scalar(theta,phi) result(er)
    real(wp), intent(in) :: theta,phi
    real(wp), dimension(3) :: er
  
    ! fixme: error checking
    er(1)=sin(theta)*cos(phi)
    er(2)=sin(theta)*sin(phi)
    er(3)=cos(theta)
  end function er_spherical_scalar
  
  !> rank 3 array:  radial unit vector (expressed in ECEF cartesian coodinates, components permuted as ix,iy,iz)
  pure function er_spherical_rank3(theta,phi) result(er)
    real(wp), dimension(:,:,:), intent(in) :: theta,phi
    real(wp), dimension(1:size(theta,1),1:size(theta,2),1:size(theta,3),3) :: er
  
    ! fixme: error checking
    er(:,:,:,1)=sin(theta)*cos(phi)
    er(:,:,:,2)=sin(theta)*sin(phi)
    er(:,:,:,3)=cos(theta)
  end function er_spherical_rank3
  
  !> scalar:  zenith angle unit vector (expressed in ECEF cartesian coodinates
  pure function etheta_spherical_scalar(theta,phi) result(etheta)
    real(wp), intent(in) :: theta,phi
    real(wp), dimension(3) :: etheta
  
    ! fixme: error checking
    etheta(1)=cos(theta)*cos(phi)
    etheta(2)=cos(theta)*sin(phi)
    etheta(3)=-sin(theta)
  end function etheta_spherical_scalar
  
  !> rank3:  zenith angle unit vector (expressed in ECEF cartesian coodinates
  pure function etheta_spherical_rank3(theta,phi) result(etheta)
    real(wp), dimension(:,:,:), intent(in) :: theta,phi
    real(wp), dimension(1:size(theta,1),1:size(theta,2),1:size(theta,3),3) :: etheta
  
    ! fixme: error checking
    etheta(:,:,:,1)=cos(theta)*cos(phi)
    etheta(:,:,:,2)=cos(theta)*sin(phi)
    etheta(:,:,:,3)=-sin(theta)
  end function etheta_spherical_rank3

  !> sclara:  azimuth angle unit vector (ECEF cart.)
  pure function ephi_spherical_scalar(theta,phi) result(ephi)
    real(wp), intent(in) :: theta,phi
    real(wp), dimension(3) :: ephi
  
    ! fixme: error checking
    ephi(1)=-sin(phi)
    ephi(2)=cos(phi)
    ephi(3)=0
  end function ephi_spherical_scalar
  
  !> rank3:  azimuth angle unit vector (ECEF cart.)
  pure function ephi_spherical_rank3(theta,phi) result(ephi)
    real(wp), dimension(:,:,:), intent(in) :: theta,phi
    real(wp), dimension(1:size(theta,1),1:size(theta,2),1:size(theta,3),3) :: ephi
  
    ! fixme: error checking
    ephi(:,:,:,1)=-sin(phi)
    ephi(:,:,:,2)=cos(phi)
    ephi(:,:,:,3)=0
  end function ephi_spherical_rank3
end module spherical
