submodule(meshobj_dipole) dipole_fns

!> This submodule contains the functions for which we need to find roots in order to transform
!   dipole to spherical coordinates

implicit none

contains

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
