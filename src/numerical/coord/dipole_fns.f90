submodule(meshobj_dipole) dipole_fns

!> This submodule contains the functions for which we need to find roots in order to transform
!   dipole to spherical coordinates.  This should arguably be its own standalone module...

implicit none

! magnetic pole location in geographic coordinates
real(wp), parameter :: thetan=11*pi/180
real(wp), parameter :: phin=289*pi/180

! options structure for Newton iterations
type(newtopts) :: newtparms

contains

!> convert a single q,p pair into r,theta
module procedure qp2rtheta
  real(wp), dimension(2) :: parms
  real(wp) :: r0
  procedure(objfun), pointer :: f
  procedure(objfun_deriv), pointer :: fprime
  integer :: maxrestart, maxr, r0step
  integer :: it,ir0
  logical :: converged

  ! Set parameters of the restart and Newton iterations
  maxrestart=400
  maxr=100*Re
  r0step=0.25*Re
  newtparms%maxit=100
  newtparms%derivtol=1e-18
  newtparms%tol=1e-11
  newtparms%verbose=.false.
  f=>rpoly
  fprime=>rpoly_deriv
  parms=[q,p]

  ! Newton iterations with restarting (see parameters above for limits) until we get a satisfactory result
  r=0; converged=.false.; ir0=1;
  do while (.not. converged .and. ir0<maxrestart .and. (r<=0 .or. r>maxr))
    r0=(ir0-1)*(r0step)    ! change starting point in increments of 0.25 Re until we get a "good" answer
    call newton_exact(f,fprime,r0,parms,newtparms,r,it,converged)
    ir0=ir0+1
  end do

  ! Once we have r can algebraically solve for theta
  theta=qr2theta(q,r)
end procedure qp2rtheta


!> convert a single point r,theta to q,p
module procedure rtheta2qp
  q=Re**2/r**2*cos(theta)
  p=r/Re/(sin(theta)**2)
end procedure rtheta2qp


!> find theta given q,r
module procedure qr2theta
  theta=acos(q*(r/Re)**2)
end procedure qr2theta


!> convert geographic coordinates to geomagnetic; do not use at the magnetic pole!!!
module procedure geog2geomag_rank3
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
end procedure geog2geomag_rank3
module procedure geog2geomag_scalar
  real(wp) :: glonwrap
  real(wp) :: thetag
  real(wp) :: phig
  real(wp) :: argtmp,alpha

  glonwrap=mod(glon,360._wp)
  thetag=pi/2._wp-glat*pi/180._wp
  phig=glonwrap*pi/180._wp

  theta = acos(cos(thetag)*cos(thetan)+sin(thetag)*sin(thetan)*cos(phig-phin))
  argtmp = (cos(thetag)-cos(theta)*cos(thetan))/(sin(theta)*sin(thetan))
  alpha = acos( max(min(argtmp,1._wp),-1._wp) )

  if (phin>phig .and. phin-phig>pi .or. phin<phig .and. phig-phin<pi) then
    phi=pi-alpha
  else
    phi=alpha+pi
  end if
end procedure geog2geomag_scalar


!> convert geomagnetic coordinates to geographic
module procedure geomag2geog_rank3
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
end procedure geomag2geog_rank3
module procedure geomag2geog_scalar
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
end procedure geomag2geog_scalar


!> convert geocentric distance into altitude (assume spherical Earth but could use other model)
module procedure r2alt
  alt=r-Re
end procedure r2alt


!> convert altitude to geocentric distance
module procedure alt2r
  r=alt+Re
end procedure alt2r


!> objective function for newton iterations for solutions of roots for r
module procedure rpoly
  real(wp) ::  q,p

  q=parms(1); p=parms(2);
  fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
end procedure rpoly


!> derivative objective function for newton iterations for roots of r
module procedure rpoly_deriv
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
end procedure rpoly_deriv

end submodule dipole_fns
