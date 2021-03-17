module dipole

!> Contains data and subroutines for converting between dipole coordinates and other coordinate systems

use phys_consts, only : Re,pi

real(wp), private, parameter :: thetan=11*pi/180
real(wp), private, parameter :: phin=289*pi/180

contains

subroutine qp2rtheta(q,p,r,theta)

end subroutine qp2rtheta

subroutine rtheta2qp(r,theta,q,p)

end subroutine rtheta2qp

function fval=objfunr(r,parms)

end function objfunr

function fpval=objfunr_deriv(r,parms)


end module dipole
