submodule(meshobj_dipole) dipole_fns

!> This submodule contains the functions for which we need to find roots in order to transform
!   dipole to spherical coordinates

implicit none

contains

!> objective function for newton iterations for solutions of roots for r
!real(wp) function rpoly(x,parms) result(fval)
!  real(wp), intent(in) :: x
!  real(wp), dimension(:), intent(in) :: parms
module procedure rpoly
  real(wp) ::  q,p

  q=parms(1); p=parms(2);
  fval=q**2*(x/Re)**4 + 1/p*(x/Re) - 1
!end function rpoly
end procedure rpoly


!> derivative objective function for newton iterations for roots of r
!real(wp) function rpoly_deriv(x,parms) result(fval_deriv)
!  real(wp), intent(in) :: x
!  real(wp), dimension(:), intent(in) :: parms
module procedure rpoly_deriv
  real(wp) :: q,p

  q=parms(1); p=parms(2);
  fval_deriv=4/Re*q**2*(x/Re)**3 + 1/p/Re
!end function rpoly_deriv
end procedure rpoly_deriv

end submodule dipole_fns
