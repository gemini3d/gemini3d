module newton

!> uses
use phys_consts, only : wp

!> a structure def. for containing the options for the newton method procedure
type newtopts
  real(wp) :: derivtol=1e-18
  integer :: maxit=100
  real(wp) :: tol=1e-9
  logical :: verbose=.false.
end type newtopts

!> these interfaced define the types of functions needed to run Newton iterations
!   need to match custom functions defined in using modules
abstract interface
  function objfun(x,parms)
    real(wp), intent(in) :: x
    real(wp), dimension(:), intent(in) :: parms
  end function objfun
end abstract interface

abstract interface
  function objfun_deriv(x,parms)
    real(wp), intent(:) :: x
    real(wp), dimension(:), intent(in) :: parms
  end function objfun_deriv
end abstract interface

contains

!> this implmements the exact Newton method for solving a nonlinear equation
subroutine newton_exact(f,fprime,x0,parms,newtparms,root,it,converged)
  procedure(objfun), pointer :: f
  procedure(objfun_deriv), pointer :: fprime
  real(wp) :: x0                      ! starting point for newton iteration
  real(wp),dimension(:) :: parms      ! fixed parameters of the newton iteration, f,fprime must accommodate whatever size array is passed in
  type(newtopts) :: newtparms         ! options for the iteration that can be set by the user
  integer :: it=0
  logical :: converged=.false.
  real(wp) :: derivative


  ! check starting point is not too close to inflection
  if (abs(fprime(x0,parms))<newtparm%dervitol) then
    print*, 'Warning:  starting near inflection point, please change initial guess!'
    return
  end if

  ! Newton iteration main loop
  it=1; root=c0; fval=f(root,parms);
  do while (.not. converged .and. it <= newtparms%maxit)
    derivative=fprime(root,parms)
    if (abs(derivative)<newtparm%dervitol) then
      print*, 'Warning:  Encountered inflection point during iteration:  ',it
      return
    else
      root=root-fval/derivative
      fval=f(root,parms)
      if (verbose) then
        print*, ' Iteration ',it,'; root ',root,' fval ',fval,' derivative ',derivative
      end if
      it=it+1
      converged=abs(fval)<newtparms%tol
    end if
  end do
  it=it-1
  return

end subroutine newton_exact

end module newton
