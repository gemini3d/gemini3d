module coordinate_transforms

!! Analytical spherical-Earth coordinate transformations.
!!
!! Angular arguments are radians and linear arguments are meters. Magnetic
!! dipole q and p use GEMINI's dimensionless normalization:
!!
!!   q = (Re/r)**2 * cos(theta)
!!   p = (r/Re) / sin(theta)**2
!!
!! This module is intentionally independent of GEMINI's existing production
!! coordinate paths so it can be validated before any caller is migrated.

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
use phys_consts, only : wp, Re, pi
use geomagnetic, only : geog2geomag, geomag2geog

implicit none (type, external)

private
public :: geographic_to_ecef, ecef_to_geographic
public :: ecef_to_enu, enu_to_ecef
public :: enu_to_aer, aer_to_enu, ecef_to_aer, aer_to_ecef
public :: geographic_to_geomagnetic, geomagnetic_to_geographic
public :: spherical_to_magdip, magdip_to_spherical
public :: geocentric_to_magdip, magdip_to_geocentric
public :: geomagnetic_ecef_to_magdip, magdip_to_geomagnetic_ecef
public :: geographic_ecef_to_magdip, magdip_to_geographic_ecef
public :: angle_difference

real(wp), parameter :: two_pi = 2._wp*pi

contains

elemental subroutine geographic_to_ecef(lat, lon, alt, x, y, z)
  !! Spherical geographic latitude, longitude, altitude to geographic ECEF.
  real(wp), intent(in) :: lat, lon, alt
  real(wp), intent(out) :: x, y, z
  real(wp) :: radius

  radius = Re + alt
  x = radius*cos(lat)*cos(lon)
  y = radius*cos(lat)*sin(lon)
  z = radius*sin(lat)
end subroutine geographic_to_ecef


elemental subroutine ecef_to_geographic(x, y, z, lat, lon, alt)
  !! Geographic ECEF to spherical geographic latitude, longitude, altitude.
  real(wp), intent(in) :: x, y, z
  real(wp), intent(out) :: lat, lon, alt
  real(wp) :: horizontal, radius, nan

  radius = sqrt(x*x + y*y + z*z)
  if (radius == 0._wp) then
    nan = quiet_nan()
    lat = nan
    lon = nan
    alt = nan
    return
  end if

  horizontal = sqrt(x*x + y*y)
  lat = atan2(z, horizontal)
  if (horizontal == 0._wp) then
    lon = 0._wp
  else
    lon = wrap_to_pi(atan2(y, x))
  end if
  alt = radius - Re
end subroutine ecef_to_geographic


elemental subroutine ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt, east, north, up)
  !! Geographic ECEF to local east, north, up at a spherical reference point.
  real(wp), intent(in) :: x, y, z
  real(wp), intent(in) :: ref_lat, ref_lon, ref_alt
  real(wp), intent(out) :: east, north, up
  real(wp) :: ref_x, ref_y, ref_z, dx, dy, dz

  call geographic_to_ecef(ref_lat, ref_lon, ref_alt, ref_x, ref_y, ref_z)
  dx = x - ref_x
  dy = y - ref_y
  dz = z - ref_z

  east = -sin(ref_lon)*dx + cos(ref_lon)*dy
  north = -sin(ref_lat)*cos(ref_lon)*dx - sin(ref_lat)*sin(ref_lon)*dy + cos(ref_lat)*dz
  up = cos(ref_lat)*cos(ref_lon)*dx + cos(ref_lat)*sin(ref_lon)*dy + sin(ref_lat)*dz
end subroutine ecef_to_enu


elemental subroutine enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_alt, x, y, z)
  !! Local east, north, up to geographic ECEF at a spherical reference point.
  real(wp), intent(in) :: east, north, up
  real(wp), intent(in) :: ref_lat, ref_lon, ref_alt
  real(wp), intent(out) :: x, y, z
  real(wp) :: ref_x, ref_y, ref_z, dx, dy, dz

  call geographic_to_ecef(ref_lat, ref_lon, ref_alt, ref_x, ref_y, ref_z)

  dx = -sin(ref_lon)*east - sin(ref_lat)*cos(ref_lon)*north + cos(ref_lat)*cos(ref_lon)*up
  dy = cos(ref_lon)*east - sin(ref_lat)*sin(ref_lon)*north + cos(ref_lat)*sin(ref_lon)*up
  dz = cos(ref_lat)*north + sin(ref_lat)*up

  x = ref_x + dx
  y = ref_y + dy
  z = ref_z + dz
end subroutine enu_to_ecef


elemental subroutine enu_to_aer(east, north, up, azimuth, elevation, slant_range)
  !! ENU to azimuth clockwise from north, elevation, and slant range.
  real(wp), intent(in) :: east, north, up
  real(wp), intent(out) :: azimuth, elevation, slant_range
  real(wp) :: horizontal

  horizontal = sqrt(east*east + north*north)
  slant_range = sqrt(horizontal*horizontal + up*up)
  if (slant_range == 0._wp) then
    azimuth = 0._wp
    elevation = 0._wp
    return
  end if

  azimuth = modulo(atan2(east, north), two_pi)
  elevation = atan2(up, horizontal)
end subroutine enu_to_aer


elemental subroutine aer_to_enu(azimuth, elevation, slant_range, east, north, up)
  !! Azimuth clockwise from north, elevation, and slant range to ENU.
  real(wp), intent(in) :: azimuth, elevation, slant_range
  real(wp), intent(out) :: east, north, up
  real(wp) :: nan

  if (slant_range < 0._wp .or. .not. ieee_is_finite(slant_range)) then
    nan = quiet_nan()
    east = nan
    north = nan
    up = nan
    return
  end if

  east = slant_range*cos(elevation)*sin(azimuth)
  north = slant_range*cos(elevation)*cos(azimuth)
  up = slant_range*sin(elevation)
end subroutine aer_to_enu


elemental subroutine ecef_to_aer(x, y, z, ref_lat, ref_lon, ref_alt, azimuth, elevation, slant_range)
  !! Geographic ECEF to AER through the reciprocal local ENU transform.
  real(wp), intent(in) :: x, y, z
  real(wp), intent(in) :: ref_lat, ref_lon, ref_alt
  real(wp), intent(out) :: azimuth, elevation, slant_range
  real(wp) :: east, north, up

  call ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt, east, north, up)
  call enu_to_aer(east, north, up, azimuth, elevation, slant_range)
end subroutine ecef_to_aer


elemental subroutine aer_to_ecef(azimuth, elevation, slant_range, ref_lat, ref_lon, ref_alt, x, y, z)
  !! AER to geographic ECEF through the reciprocal local ENU transform.
  real(wp), intent(in) :: azimuth, elevation, slant_range
  real(wp), intent(in) :: ref_lat, ref_lon, ref_alt
  real(wp), intent(out) :: x, y, z
  real(wp) :: east, north, up

  call aer_to_enu(azimuth, elevation, slant_range, east, north, up)
  call enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_alt, x, y, z)
end subroutine aer_to_ecef


elemental subroutine geographic_to_geomagnetic(lat, lon, theta, phi)
  !! Geographic latitude/longitude to GEMINI magnetic colatitude/longitude.
  real(wp), intent(in) :: lat, lon
  real(wp), intent(out) :: theta, phi

  ! geog2geomag uses MOD internally, so normalize negative longitudes first.
  call geog2geomag(modulo(lon, two_pi)*180._wp/pi, lat*180._wp/pi, phi, theta)
  phi = modulo(phi, two_pi)
end subroutine geographic_to_geomagnetic


elemental subroutine geomagnetic_to_geographic(theta, phi, lat, lon)
  !! GEMINI magnetic colatitude/longitude to geographic latitude/longitude.
  real(wp), intent(in) :: theta, phi
  real(wp), intent(out) :: lat, lon
  real(wp) :: glat_deg, glon_deg

  ! geomag2geog also uses MOD internally; accept either longitude convention.
  call geomag2geog(modulo(phi, two_pi), theta, glon_deg, glat_deg)
  lat = glat_deg*pi/180._wp
  lon = wrap_to_pi(glon_deg*pi/180._wp)
end subroutine geomagnetic_to_geographic


elemental subroutine spherical_to_magdip(r, theta, phi, q, p, phi_out)
  !! Magnetic spherical radius, colatitude, longitude to normalized q, p, phi.
  real(wp), intent(in) :: r, theta, phi
  real(wp), intent(out) :: q, p, phi_out
  real(wp) :: sin_theta, nan

  if (r <= 0._wp .or. .not. ieee_is_finite(r)) then
    nan = quiet_nan()
    q = nan
    p = nan
    phi_out = nan
    return
  end if

  sin_theta = sin(theta)
  q = (Re/r)**2*cos(theta)
  p = (r/Re)/(sin_theta*sin_theta)
  phi_out = wrap_to_pi(phi)
end subroutine spherical_to_magdip


elemental subroutine magdip_to_spherical(q, p, phi, r, theta, phi_out)
  !! Analytical normalized dipole inverse following Swisdak (2006).
  real(wp), intent(in) :: q, p, phi
  real(wp), intent(out) :: r, theta, phi_out
  real(wp) :: alpha, beta, gamma, mu, r_over_re, sin_theta_sq
  real(wp) :: theta_candidate, nan

  if (p <= 0._wp .or. .not. ieee_is_finite(q) .or. .not. ieee_is_finite(p)) then
    nan = quiet_nan()
    r = nan
    theta = nan
    phi_out = nan
    return
  end if

  alpha = (256._wp/27._wp)*q*q*p**4
  beta = (1._wp + sqrt(1._wp + alpha))**(2._wp/3._wp)
  gamma = alpha**(1._wp/3._wp)
  mu = 0.5_wp*((beta*beta + beta*gamma + gamma*gamma)/beta)**1.5_wp

  r_over_re = 4._wp*mu*p / ((1._wp + mu)*(1._wp + sqrt(max(2._wp*mu - 1._wp, 0._wp))))
  r = Re*r_over_re
  sin_theta_sq = min(max(r_over_re/p, 0._wp), 1._wp)
  theta_candidate = asin(sqrt(sin_theta_sq))
  theta = merge(pi - theta_candidate, theta_candidate, q < 0._wp)
  phi_out = wrap_to_pi(phi)
end subroutine magdip_to_spherical


elemental subroutine geocentric_to_magdip(r, lat, lon, q, p, phi)
  !! Magnetic geocentric radius, latitude, longitude to normalized dipole.
  real(wp), intent(in) :: r, lat, lon
  real(wp), intent(out) :: q, p, phi

  call spherical_to_magdip(r, pi/2._wp - lat, lon, q, p, phi)
end subroutine geocentric_to_magdip


elemental subroutine magdip_to_geocentric(q, p, phi, r, lat, lon)
  !! Normalized dipole to magnetic geocentric radius, latitude, longitude.
  real(wp), intent(in) :: q, p, phi
  real(wp), intent(out) :: r, lat, lon
  real(wp) :: theta

  call magdip_to_spherical(q, p, phi, r, theta, lon)
  lat = pi/2._wp - theta
end subroutine magdip_to_geocentric


elemental subroutine geomagnetic_ecef_to_magdip(x, y, z, q, p, phi)
  !! Dipole-axis-aligned ECEF Cartesian coordinates to normalized dipole.
  real(wp), intent(in) :: x, y, z
  real(wp), intent(out) :: q, p, phi
  real(wp) :: r, theta, nan

  r = sqrt(x*x + y*y + z*z)
  if (r == 0._wp) then
    nan = quiet_nan()
    q = nan
    p = nan
    phi = nan
    return
  end if

  theta = atan2(sqrt(x*x + y*y), z)
  call spherical_to_magdip(r, theta, atan2(y, x), q, p, phi)
end subroutine geomagnetic_ecef_to_magdip


elemental subroutine magdip_to_geomagnetic_ecef(q, p, phi, x, y, z)
  !! Normalized dipole to dipole-axis-aligned ECEF Cartesian coordinates.
  real(wp), intent(in) :: q, p, phi
  real(wp), intent(out) :: x, y, z
  real(wp) :: r, theta, phi_sph

  call magdip_to_spherical(q, p, phi, r, theta, phi_sph)
  x = r*sin(theta)*cos(phi_sph)
  y = r*sin(theta)*sin(phi_sph)
  z = r*cos(theta)
end subroutine magdip_to_geomagnetic_ecef


elemental subroutine geographic_ecef_to_magdip(x, y, z, q, p, phi)
  !! Geographic ECEF Cartesian coordinates to GEMINI normalized dipole.
  real(wp), intent(in) :: x, y, z
  real(wp), intent(out) :: q, p, phi
  real(wp) :: lat, lon, alt, theta, phi_mag

  call ecef_to_geographic(x, y, z, lat, lon, alt)
  call geographic_to_geomagnetic(lat, lon, theta, phi_mag)
  call spherical_to_magdip(Re + alt, theta, phi_mag, q, p, phi)
end subroutine geographic_ecef_to_magdip


elemental subroutine magdip_to_geographic_ecef(q, p, phi, x, y, z)
  !! GEMINI normalized dipole to geographic ECEF Cartesian coordinates.
  real(wp), intent(in) :: q, p, phi
  real(wp), intent(out) :: x, y, z
  real(wp) :: r, theta, phi_sph, lat, lon

  call magdip_to_spherical(q, p, phi, r, theta, phi_sph)
  call geomagnetic_to_geographic(theta, phi_sph, lat, lon)
  call geographic_to_ecef(lat, lon, r - Re, x, y, z)
end subroutine magdip_to_geographic_ecef


elemental function angle_difference(a, b) result(difference)
  !! Minimal signed angular difference a-b in [-pi, pi).
  real(wp), intent(in) :: a, b
  real(wp) :: difference

  difference = wrap_to_pi(a - b)
end function angle_difference


elemental function wrap_to_pi(angle) result(wrapped)
  real(wp), intent(in) :: angle
  real(wp) :: wrapped

  wrapped = modulo(angle + pi, two_pi) - pi
end function wrap_to_pi


elemental function quiet_nan() result(nan)
  real(wp) :: nan

  nan = ieee_value(0._wp, ieee_quiet_nan)
end function quiet_nan

end module coordinate_transforms
