program coordinate_transforms_testdriver

use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use phys_consts, only : wp, Re, pi
use geomagnetic, only : set_magnetic_pole
use dipole, only : qp2rtheta, qp2rtheta_newton
use coordinate_transforms, only : geographic_to_ecef, ecef_to_geographic, &
  ecef_to_enu, enu_to_ecef, enu_to_aer, aer_to_enu, ecef_to_aer, aer_to_ecef, &
  geographic_to_geomagnetic, geomagnetic_to_geographic, &
  spherical_to_magdip, magdip_to_spherical, geocentric_to_magdip, magdip_to_geocentric, &
  geomagnetic_ecef_to_magdip, magdip_to_geomagnetic_ecef, &
  geographic_ecef_to_magdip, magdip_to_geographic_ecef, angle_difference

implicit none (type, external)

real(wp), parameter :: linear_tol = 1.e-6_wp
real(wp), parameter :: angle_tol = 1.e-12_wp
real(wp), parameter :: dipole_angle_tol = 1.e-9_wp

call test_geographic_ecef_controls()
call test_enu_aer_controls()
call test_scalar_roundtrips()
call test_array_roundtrips()
call test_geomagnetic_roundtrips()
call test_dipole_roundtrips()
call test_invalid_inputs()
call compare_analytic_with_newton()

print '(a)', 'Analytical coordinate transformation tests: PASS'

contains

subroutine test_geographic_ecef_controls()
  real(wp) :: x, y, z, lat, lon, alt

  call geographic_to_ecef(0._wp, 0._wp, 0._wp, x, y, z)
  call assert_close('equator x', x, Re, linear_tol)
  call assert_close('equator y', y, 0._wp, linear_tol)
  call assert_close('equator z', z, 0._wp, linear_tol)

  call geographic_to_ecef(0._wp, pi/2._wp, 0._wp, x, y, z)
  call assert_close('east equator x', x, 0._wp, linear_tol)
  call assert_close('east equator y', y, Re, linear_tol)

  call geographic_to_ecef(pi/2._wp, 0._wp, 0._wp, x, y, z)
  call assert_close('north pole x', x, 0._wp, linear_tol)
  call assert_close('north pole z', z, Re, linear_tol)

  call geographic_to_ecef(0._wp, pi, 0._wp, x, y, z)
  call assert_close('dateline x', x, -Re, linear_tol)
  call ecef_to_geographic(x, y, z, lat, lon, alt)
  call assert_close('dateline latitude', lat, 0._wp, angle_tol)
  call assert_angle_close('dateline longitude', lon, -pi, angle_tol)
  call assert_close('dateline altitude', alt, 0._wp, linear_tol)

  call ecef_to_geographic(0._wp, 0._wp, Re + 1000._wp, lat, lon, alt)
  call assert_close('pole inverse latitude', lat, pi/2._wp, angle_tol)
  call assert_close('pole inverse longitude convention', lon, 0._wp, angle_tol)
  call assert_close('pole inverse altitude', alt, 1000._wp, linear_tol)
end subroutine test_geographic_ecef_controls


subroutine test_enu_aer_controls()
  real(wp) :: east, north, up, azimuth, elevation, slant_range
  real(wp), parameter :: hundred = 100._wp

  call ecef_to_enu(Re, hundred, 0._wp, 0._wp, 0._wp, 0._wp, east, north, up)
  call assert_close('ENU east axis', east, hundred, linear_tol)
  call assert_close('ENU east north', north, 0._wp, linear_tol)
  call assert_close('ENU east up', up, 0._wp, linear_tol)

  call ecef_to_enu(Re, 0._wp, hundred, 0._wp, 0._wp, 0._wp, east, north, up)
  call assert_close('ENU north axis', north, hundred, linear_tol)

  call ecef_to_enu(Re + hundred, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, east, north, up)
  call assert_close('ENU up axis', up, hundred, linear_tol)

  call enu_to_aer(hundred, 0._wp, 0._wp, azimuth, elevation, slant_range)
  call assert_close('AER east azimuth', azimuth, pi/2._wp, angle_tol)
  call assert_close('AER east elevation', elevation, 0._wp, angle_tol)
  call assert_close('AER east range', slant_range, hundred, linear_tol)

  call enu_to_aer(0._wp, hundred, 0._wp, azimuth, elevation, slant_range)
  call assert_close('AER north azimuth', azimuth, 0._wp, angle_tol)

  call enu_to_aer(0._wp, 0._wp, hundred, azimuth, elevation, slant_range)
  call assert_close('AER up elevation', elevation, pi/2._wp, angle_tol)

  call aer_to_enu(pi/2._wp, 0._wp, hundred, east, north, up)
  call assert_close('AER inverse east', east, hundred, linear_tol)
  call assert_close('AER inverse north', north, 0._wp, linear_tol)
end subroutine test_enu_aer_controls


subroutine test_scalar_roundtrips()
  real(wp) :: lat, lon, alt, lat_back, lon_back, alt_back
  real(wp) :: ref_lat, ref_lon, ref_alt, x, y, z
  real(wp) :: east, north, up, east_back, north_back, up_back
  real(wp) :: azimuth, elevation, slant_range, azimuth_back, elevation_back, range_back
  real(wp) :: theta, phi, theta_back, phi_back, phi_dipole, radius, radius_back, q, p

  lat = 0.37_wp
  lon = -1.24_wp
  alt = 123456._wp
  call geographic_to_ecef(lat, lon, alt, x, y, z)
  call ecef_to_geographic(x, y, z, lat_back, lon_back, alt_back)
  call assert_close('scalar geographic latitude', lat_back, lat, angle_tol)
  call assert_angle_close('scalar geographic longitude', lon_back, lon, angle_tol)
  call assert_close('scalar geographic altitude', alt_back, alt, linear_tol)

  ref_lat = 0.58_wp
  ref_lon = -1.70_wp
  ref_alt = 1250._wp
  east = 12345._wp
  north = -67890._wp
  up = 456._wp
  call enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_alt, x, y, z)
  call ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt, east_back, north_back, up_back)
  call assert_close('scalar ECEF/ENU east', east_back, east, linear_tol)
  call assert_close('scalar ECEF/ENU north', north_back, north, linear_tol)
  call assert_close('scalar ECEF/ENU up', up_back, up, linear_tol)

  call enu_to_aer(east, north, up, azimuth, elevation, slant_range)
  call aer_to_enu(azimuth, elevation, slant_range, east_back, north_back, up_back)
  call assert_close('scalar ENU/AER east', east_back, east, linear_tol)
  call assert_close('scalar ENU/AER north', north_back, north, linear_tol)
  call assert_close('scalar ENU/AER up', up_back, up, linear_tol)
  call aer_to_ecef(azimuth, elevation, slant_range, ref_lat, ref_lon, ref_alt, x, y, z)
  call ecef_to_aer(x, y, z, ref_lat, ref_lon, ref_alt, azimuth_back, elevation_back, range_back)
  call assert_angle_close('scalar ECEF/AER azimuth', azimuth_back, azimuth, angle_tol)
  call assert_close('scalar ECEF/AER elevation', elevation_back, elevation, angle_tol)
  call assert_close('scalar ECEF/AER range', range_back, slant_range, linear_tol)

  call set_magnetic_pole(2024)
  call geographic_to_geomagnetic(lat, lon, theta, phi)
  call geomagnetic_to_geographic(theta, phi, lat_back, lon_back)
  call assert_close('scalar geomagnetic latitude', lat_back, lat, angle_tol)
  call assert_angle_close('scalar geomagnetic longitude', lon_back, lon, angle_tol)

  radius = Re + alt
  call spherical_to_magdip(radius, theta, phi, q, p, phi_dipole)
  call magdip_to_spherical(q, p, phi_dipole, radius_back, theta_back, phi_back)
  call assert_close('scalar spherical/dipole radius', radius_back, radius, linear_tol)
  call assert_close('scalar spherical/dipole theta', theta_back, theta, angle_tol)
  call assert_angle_close('scalar spherical/dipole phi', phi_back, phi, angle_tol)

  call geocentric_to_magdip(radius, pi/2._wp-theta, phi, q, p, phi_back)
  call magdip_to_geocentric(q, p, phi_back, radius_back, lat_back, lon_back)
  call assert_close('scalar geocentric/dipole radius', radius_back, radius, linear_tol)
  call assert_close('scalar geocentric/dipole latitude', lat_back, pi/2._wp-theta, angle_tol)
  call assert_angle_close('scalar geocentric/dipole longitude', lon_back, phi, angle_tol)

  x = radius*sin(theta)*cos(phi)
  y = radius*sin(theta)*sin(phi)
  z = radius*cos(theta)
  call geomagnetic_ecef_to_magdip(x, y, z, q, p, phi_back)
  call magdip_to_geomagnetic_ecef(q, p, phi_back, east, north, up)
  call assert_close('scalar magnetic ECEF/dipole x', east, x, linear_tol)
  call assert_close('scalar magnetic ECEF/dipole y', north, y, linear_tol)
  call assert_close('scalar magnetic ECEF/dipole z', up, z, linear_tol)

  call geographic_to_ecef(lat, lon, alt, x, y, z)
  call geographic_ecef_to_magdip(x, y, z, q, p, phi_back)
  call magdip_to_geographic_ecef(q, p, phi_back, east, north, up)
  call assert_close('scalar geographic ECEF/dipole x', east, x, linear_tol)
  call assert_close('scalar geographic ECEF/dipole y', north, y, linear_tol)
  call assert_close('scalar geographic ECEF/dipole z', up, z, linear_tol)
  call set_magnetic_pole(0)
end subroutine test_scalar_roundtrips


subroutine test_array_roundtrips()
  integer, parameter :: nlat=7, nlon=9, nalt=4
  integer :: i, j, k
  real(wp) :: lat(nlat,nlon,nalt), lon(nlat,nlon,nalt), alt(nlat,nlon,nalt)
  real(wp) :: lat_back(nlat,nlon,nalt), lon_back(nlat,nlon,nalt), alt_back(nlat,nlon,nalt)
  real(wp) :: x(nlat,nlon,nalt), y(nlat,nlon,nalt), z(nlat,nlon,nalt)
  real(wp) :: east(nlat,nlon,nalt), north(nlat,nlon,nalt), up(nlat,nlon,nalt)
  real(wp) :: east_back(nlat,nlon,nalt), north_back(nlat,nlon,nalt), up_back(nlat,nlon,nalt)
  real(wp) :: azimuth(nlat,nlon,nalt), elevation(nlat,nlon,nalt), slant_range(nlat,nlon,nalt)
  real(wp) :: ref_lat, ref_lon, ref_alt

  do k=1,nalt
    do j=1,nlon
      do i=1,nlat
        lat(i,j,k) = (-80._wp + real(i-1,wp)*160._wp/real(nlat-1,wp))*pi/180._wp
        lon(i,j,k) = (-170._wp + real(j-1,wp)*340._wp/real(nlon-1,wp))*pi/180._wp
        alt(i,j,k) = -1000._wp + real(k-1,wp)*501000._wp/real(nalt-1,wp)
      end do
    end do
  end do

  call geographic_to_ecef(lat, lon, alt, x, y, z)
  call ecef_to_geographic(x, y, z, lat_back, lon_back, alt_back)
  call assert_max('geographic latitude array', abs(lat_back-lat), angle_tol)
  call assert_max('geographic longitude array', abs(angle_difference(lon_back,lon)), angle_tol)
  call assert_max('geographic altitude array', abs(alt_back-alt), linear_tol)

  ref_lat = 33.3444_wp*pi/180._wp
  ref_lon = -97.5702_wp*pi/180._wp
  ref_alt = 250._wp
  do k=1,nalt
    do j=1,nlon
      do i=1,nlat
        east(i,j,k) = -500000._wp + real(i-1,wp)*1000000._wp/real(nlat-1,wp)
        north(i,j,k) = -700000._wp + real(j-1,wp)*1400000._wp/real(nlon-1,wp)
        up(i,j,k) = -1000._wp + real(k-1,wp)*501000._wp/real(nalt-1,wp)
      end do
    end do
  end do

  call enu_to_ecef(east, north, up, ref_lat, ref_lon, ref_alt, x, y, z)
  call ecef_to_enu(x, y, z, ref_lat, ref_lon, ref_alt, east_back, north_back, up_back)
  call assert_max('ENU east array', abs(east_back-east), linear_tol)
  call assert_max('ENU north array', abs(north_back-north), linear_tol)
  call assert_max('ENU up array', abs(up_back-up), linear_tol)

  call enu_to_aer(east, north, up, azimuth, elevation, slant_range)
  call aer_to_enu(azimuth, elevation, slant_range, east_back, north_back, up_back)
  call assert_max('AER east array', abs(east_back-east), linear_tol)
  call assert_max('AER north array', abs(north_back-north), linear_tol)
  call assert_max('AER up array', abs(up_back-up), linear_tol)

  call aer_to_ecef(azimuth, elevation, slant_range, ref_lat, ref_lon, ref_alt, x, y, z)
  call ecef_to_aer(x, y, z, ref_lat, ref_lon, ref_alt, azimuth, elevation, slant_range)
  call aer_to_enu(azimuth, elevation, slant_range, east_back, north_back, up_back)
  call assert_max('ECEF/AER east array', abs(east_back-east), linear_tol)
  call assert_max('ECEF/AER north array', abs(north_back-north), linear_tol)
  call assert_max('ECEF/AER up array', abs(up_back-up), linear_tol)
end subroutine test_array_roundtrips


subroutine test_geomagnetic_roundtrips()
  integer, parameter :: nlat=7, nlon=9, nalt=2
  integer :: i, j, k, year
  real(wp) :: lat(nlat,nlon,nalt), lon(nlat,nlon,nalt)
  real(wp) :: theta(nlat,nlon,nalt), phi(nlat,nlon,nalt)
  real(wp) :: lat_back(nlat,nlon,nalt), lon_back(nlat,nlon,nalt)

  do k=1,nalt
    do j=1,nlon
      do i=1,nlat
        lat(i,j,k) = (-70._wp + real(i-1,wp)*140._wp/real(nlat-1,wp))*pi/180._wp
        lon(i,j,k) = (-160._wp + real(j-1,wp)*320._wp/real(nlon-1,wp))*pi/180._wp
      end do
    end do
  end do

  do year=0,2024,2024
    call set_magnetic_pole(year)
    call geographic_to_geomagnetic(lat, lon, theta, phi)
    call geomagnetic_to_geographic(theta, phi, lat_back, lon_back)
    call assert_max('geomagnetic latitude array', abs(lat_back-lat), angle_tol)
    call assert_max('geomagnetic longitude array', abs(angle_difference(lon_back,lon)), angle_tol)
  end do
  call set_magnetic_pole(0)
end subroutine test_geomagnetic_roundtrips


subroutine test_dipole_roundtrips()
  integer, parameter :: nr=5, ntheta=9, nphi=7
  integer :: i, j, k
  real(wp) :: r(nr,ntheta,nphi), theta(nr,ntheta,nphi), phi(nr,ntheta,nphi)
  real(wp) :: q(nr,ntheta,nphi), p(nr,ntheta,nphi), phi_dip(nr,ntheta,nphi)
  real(wp) :: r_back(nr,ntheta,nphi), theta_back(nr,ntheta,nphi), phi_back(nr,ntheta,nphi)
  real(wp) :: x(nr,ntheta,nphi), y(nr,ntheta,nphi), z(nr,ntheta,nphi)
  real(wp) :: x_back(nr,ntheta,nphi), y_back(nr,ntheta,nphi), z_back(nr,ntheta,nphi)
  real(wp) :: lat(nr,ntheta,nphi), lon(nr,ntheta,nphi), alt(nr,ntheta,nphi)
  real(wp) :: lat_back(nr,ntheta,nphi)

  do k=1,nphi
    do j=1,ntheta
      do i=1,nr
        r(i,j,k) = Re*(0.9_wp + real(i-1,wp)*4.1_wp/real(nr-1,wp))
        theta(i,j,k) = 1.e-5_wp + real(j-1,wp)*(pi-2.e-5_wp)/real(ntheta-1,wp)
        phi(i,j,k) = -pi + real(k-1,wp)*two_pi()/real(nphi-1,wp)
      end do
    end do
  end do

  call spherical_to_magdip(r, theta, phi, q, p, phi_dip)
  call magdip_to_spherical(q, p, phi_dip, r_back, theta_back, phi_back)
  call assert_max('dipole radius array', abs(r_back-r), linear_tol)
  call assert_max('dipole theta interior array', abs(theta_back(:,2:ntheta-1,:)-theta(:,2:ntheta-1,:)), angle_tol)
  call assert_max('dipole theta array', abs(theta_back-theta), dipole_angle_tol)
  call assert_max('dipole phi array', abs(angle_difference(phi_back,phi)), angle_tol)

  lat = pi/2._wp - theta
  call geocentric_to_magdip(r, lat, phi, q, p, phi_dip)
  call magdip_to_geocentric(q, p, phi_dip, r_back, lat_back, phi_back)
  call assert_max('geocentric radius array', abs(r_back-r), linear_tol)
  call assert_max('geocentric latitude interior array', &
    abs(lat_back(:,2:ntheta-1,:)-lat(:,2:ntheta-1,:)), angle_tol)
  call assert_max('geocentric latitude array', abs(lat_back-lat), dipole_angle_tol)

  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(theta)
  call geomagnetic_ecef_to_magdip(x, y, z, q, p, phi_dip)
  call magdip_to_geomagnetic_ecef(q, p, phi_dip, x_back, y_back, z_back)
  call assert_max('magnetic ECEF x array', abs(x_back-x), linear_tol)
  call assert_max('magnetic ECEF y array', abs(y_back-y), linear_tol)
  call assert_max('magnetic ECEF z array', abs(z_back-z), linear_tol)

  call set_magnetic_pole(0)
  lat = (-60._wp + 120._wp*(theta-1.e-5_wp)/(pi-2.e-5_wp))*pi/180._wp
  lon = phi
  alt = r-Re
  call geographic_to_ecef(lat, lon, alt, x, y, z)
  call geographic_ecef_to_magdip(x, y, z, q, p, phi_dip)
  call magdip_to_geographic_ecef(q, p, phi_dip, x_back, y_back, z_back)
  call assert_max('geographic ECEF/dipole x array', abs(x_back-x), linear_tol)
  call assert_max('geographic ECEF/dipole y array', abs(y_back-y), linear_tol)
  call assert_max('geographic ECEF/dipole z array', abs(z_back-z), linear_tol)
end subroutine test_dipole_roundtrips


subroutine test_invalid_inputs()
  real(wp) :: a, b, c

  call ecef_to_geographic(0._wp, 0._wp, 0._wp, a, b, c)
  if (.not. all([ieee_is_nan(a),ieee_is_nan(b),ieee_is_nan(c)])) then
    error stop 'ECEF origin must return NaNs'
  end if

  call magdip_to_spherical(0._wp, 0._wp, 0._wp, a, b, c)
  if (.not. all([ieee_is_nan(a),ieee_is_nan(b),ieee_is_nan(c)])) then
    error stop 'invalid dipole p must return NaNs'
  end if

  call aer_to_enu(0._wp, 0._wp, -1._wp, a, b, c)
  if (.not. all([ieee_is_nan(a),ieee_is_nan(b),ieee_is_nan(c)])) then
    error stop 'negative AER range must return NaNs'
  end if
end subroutine test_invalid_inputs


subroutine compare_analytic_with_newton()
  integer, parameter :: nr=20, ntheta=40, nbench=100000
  integer :: i, j, n
  real(wp) :: radius, theta, q, p, phi
  real(wp) :: analytic_r, analytic_theta, analytic_phi, production_r, production_theta
  real(wp) :: newton_r, newton_theta
  real(wp) :: max_r_difference, max_theta_difference
  real(wp) :: q_back, p_back, start_time, analytic_time, newton_time
  real(wp), allocatable :: q_bench(:), p_bench(:), r_bench(:), theta_bench(:), phi_bench(:), phi_out_bench(:)

  max_r_difference = 0._wp
  max_theta_difference = 0._wp
  do j=1,ntheta
    theta = 0.05_wp + real(j-1,wp)*(pi-0.1_wp)/real(ntheta-1,wp)
    do i=1,nr
      radius = Re*(0.9_wp + real(i-1,wp)*4.1_wp/real(nr-1,wp))
      call spherical_to_magdip(radius, theta, 0._wp, q, p, phi)
      call magdip_to_spherical(q, p, phi, analytic_r, analytic_theta, analytic_phi)
      call qp2rtheta(q, p, production_r, production_theta)
      call qp2rtheta_newton(q, p, newton_r, newton_theta)
      call assert_close('production/analytic radius parity', production_r, analytic_r, linear_tol)
      call assert_close('production/analytic theta parity', production_theta, analytic_theta, angle_tol)
      max_r_difference = max(max_r_difference, abs(analytic_r-newton_r))
      max_theta_difference = max(max_theta_difference, abs(analytic_theta-newton_theta))
      call spherical_to_magdip(analytic_r, analytic_theta, 0._wp, q_back, p_back, phi)
      call assert_close('analytic q reconstruction', q_back, q, 1.e-12_wp*max(1._wp,abs(q)))
      call assert_close('analytic p reconstruction', p_back, p, 1.e-12_wp*max(1._wp,abs(p)))
    end do
  end do
  call assert_close('analytic/Newton radius parity', max_r_difference, 0._wp, 1.e-3_wp)
  call assert_close('analytic/Newton theta parity', max_theta_difference, 0._wp, 1.e-8_wp)

  allocate(q_bench(nbench), p_bench(nbench), r_bench(nbench), theta_bench(nbench), phi_bench(nbench), &
    phi_out_bench(nbench))
  do n=1,nbench
    radius = Re*(0.9_wp + 4.1_wp*real(mod(n-1,997),wp)/996._wp)
    theta = 0.05_wp + (pi-0.1_wp)*real(mod(n-1,991),wp)/990._wp
    call spherical_to_magdip(radius, theta, 0._wp, q_bench(n), p_bench(n), phi_bench(n))
  end do

  call cpu_time(start_time)
  call magdip_to_spherical(q_bench, p_bench, phi_bench, r_bench, theta_bench, phi_out_bench)
  call cpu_time(analytic_time)
  analytic_time = analytic_time-start_time

  call cpu_time(start_time)
  do n=1,nbench
    call qp2rtheta_newton(q_bench(n), p_bench(n), r_bench(n), theta_bench(n))
  end do
  call cpu_time(newton_time)
  newton_time = newton_time-start_time

  print '(a,es12.4,a)', 'Analytical dipole inverse time: ',analytic_time,' s'
  print '(a,es12.4,a)', 'Newton dipole inverse time:     ',newton_time,' s'
  if (analytic_time>0._wp) print '(a,f8.2,a)', 'Observed speed ratio:           ',newton_time/analytic_time,'x'
end subroutine compare_analytic_with_newton


subroutine assert_close(name, actual, expected, tolerance)
  character(*), intent(in) :: name
  real(wp), intent(in) :: actual, expected, tolerance

  if (abs(actual-expected)>tolerance) then
    print '(a)', trim(name)//' failed'
    print '(a,es24.16)', '  actual:    ',actual
    print '(a,es24.16)', '  expected:  ',expected
    print '(a,es24.16)', '  tolerance: ',tolerance
    error stop 'coordinate transform assertion failed'
  end if
end subroutine assert_close


subroutine assert_angle_close(name, actual, expected, tolerance)
  character(*), intent(in) :: name
  real(wp), intent(in) :: actual, expected, tolerance

  call assert_close(name, angle_difference(actual,expected), 0._wp, tolerance)
end subroutine assert_angle_close


subroutine assert_max(name, errors, tolerance)
  character(*), intent(in) :: name
  real(wp), intent(in) :: errors(:,:,:)
  real(wp), intent(in) :: tolerance

  call assert_close(name, maxval(errors), 0._wp, tolerance)
end subroutine assert_max


pure function two_pi() result(value)
  real(wp) :: value
  value = 2._wp*pi
end function two_pi

end program coordinate_transforms_testdriver
