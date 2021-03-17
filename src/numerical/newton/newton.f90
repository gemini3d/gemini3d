module newton

!> a structure for containing the options for the newton method procedure
type newtopts
  real(wp) :: derivtol=1e-18
  integer :: maxit=100
  real(wp) :: tol=1e-9
  logical :: verbose=F
end type newtopts

contains

subroutine newton_exact(f,fprime,x0,parms,newtparms)
  procedure(), pointer :: f,fprime    ! FIXME:  good form to declare procedure interface here???
  real(wp) :: x0                      ! starting point for newton iteration
  real(wp),dimension(:) :: parms      ! fixed parameters of the newton iteration, f,fprime must accommodate whatever size array is passed in
  type(newtopts) :: newtparms         ! options for the iteration that can be set by the user



end subroutine newton_exact

end module newton
