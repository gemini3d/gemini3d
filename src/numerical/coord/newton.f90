module newton

!> uses, basically we need something to tell us what precision is being used for calculations
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use phys_consts, only: wp

implicit none (type, external)
private
public :: newtopts, newton_exact, objfun, objfun_deriv

!> derived type for containing the options for the newton method procedure, by default
!   these work okay with the dipole to spherical conversion problem but can be adjusted
!   by the user for other applications
type :: newtopts
  real(wp) :: derivtol=1e-18
  integer :: maxit=100
  real(wp) :: tol=1e-9
  logical :: verbose=.false.
end type newtopts

!> these interfaced define the types of functions needed to run Newton iterations
!   need to match custom functions defined in using modules; these are defined so
!   that there are multiple objective functions (conforming to these patterns) that
!   can be used with Newton's method.  These are abstract because there are multiple
!   objective functions (in principle) that can be used with Newton's method
abstract interface
  real(wp) function objfun(x,parms) result(fval)
    import wp
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
  end function objfun
  real(wp) function objfun_deriv(x,parms) result(fval_deriv)
    import wp
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
  end function objfun_deriv
end interface


contains

!> this implmements the exact Newton method for solving a nonlinear equation
subroutine newton_exact(f,fprime,x0,parms,newtparms,root,it,converged)
  procedure(objfun), pointer, intent(in) :: f
  procedure(objfun_deriv), pointer, intent(in) :: fprime
  real(wp), intent(in) :: x0                      ! starting point for newton iteration
  real(wp),dimension(:), intent(in) :: parms      ! fixed parameters of the newton iteration, f,fprime must accommodate whatever size array is passed in
  type(newtopts), intent(in) :: newtparms         ! options for the iteration that can be set by the user
  real(wp), intent(out) :: root
  integer, intent(out) :: it
  logical, intent(out) :: converged

  real(wp) :: fval,derivative

  ! check starting point is not too close to inflection
  if (abs(fprime(x0,parms))<newtparms%derivtol) then
    write(stderr,*) 'Warning:  starting near inflection point, please change initial guess!'
    it=0; converged=.false.; root=x0;
    return
  end if

  ! Newton iteration main loop
  it=1; root=x0; fval=f(root,parms); converged=.false.
  do while (.not. converged .and. it <= newtparms%maxit)
    derivative=fprime(root,parms)
    if (abs(derivative)<newtparms%derivtol) then
      write(stderr,*) 'Warning:  Encountered inflection point during iteration:  ',it
      return
    else
      root=root-fval/derivative
      fval=f(root,parms)
      if (newtparms%verbose) then
        print*, ' Iteration ',it,'; root ',root,' fval ',fval,' derivative ',derivative
      end if
      it=it+1
      converged=abs(fval)<newtparms%tol
    end if
  end do
  it=it-1

end subroutine newton_exact

end module newton
