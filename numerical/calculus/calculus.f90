module calculus

!NEED TO MORE APPROPRIATELY NAME THE "ALT" DERIVATIVES...


!SIZES USED IN ALL DERIVATIVE PROCEDURES ARE STORED IN GRID MODULE
use phys_consts, only: wp
use grid, only : curvmesh         !we do not want the full-grid sizes (lx1,lx2,lx3) in scope since we routinely need to do subgrid derivatives
implicit none
private

!OVERLOAD ALL OF THE CALCULUS ROUTINE TO DEAL WITH THE CURVILINEAR GRID STRUCTURES
interface grad3D1
  module procedure grad3D1_curv_23
end interface grad3D1

interface grad3D2
  module procedure grad3D2_curv_23
end interface grad3D2

interface grad3D3
  module procedure grad3D3_curv_23
end interface grad3D3

interface div3D
  module procedure div3D_curv_23
end interface div3D

!ROUTINES BELOW DO NOT ACCOUNT FOR METRIC FACTORS...  As such they need to really be renamed to avoid confusion (they aren't curvilinear derivatives)

interface grad2D2
  module procedure grad2D2_curv_23
end interface grad2D2

interface grad2D3
  module procedure grad2D3_curv_23
end interface grad2D3

interface grad2D1_curv_alt
  module procedure grad2D1_curv_alt_23
end interface grad2D1_curv_alt

interface integral3D1
  module procedure integral3D1_curv
end interface integral3D1

interface integral2D1
  module procedure integral2D1_curv
end interface integral2D1

interface integral2D2
  module procedure integral2D2_curv
end interface integral2D2


public :: grad3D1, grad3D2, grad3D3, &
  grad2d1_curv_alt, grad2D2, grad2D3, grad2D3_curv_periodic, &
  div3D, &
  integral3D1, integral3D1_curv_alt, &
  chapman_a, ETD_uncoupled


contains


pure function ETD_uncoupled(f,P,L,dt)

!------------------------------------------------------------
!-------APPLY UNCOUPLED SPECIES EXPONENTIAL TIME DIFFERENCING
!-------SOURCE/LOSS ODE SOLUTION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f,P,L
real(wp), intent(in) :: dt

real(wp), dimension(size(f,1),size(f,2),size(f,3)) :: Ldt,expL

real(wp), dimension(size(f,1),size(f,2),size(f,3)) :: ETD_uncoupled

Ldt=L*dt
expL=exp(-1.0_wp*Ldt)
where (Ldt>1e-10_wp)
  ETD_uncoupled=f*expL+P/L*(1-expL)    !fast but could be inaccurate in dynamic situations
elsewhere
  ETD_uncoupled=f+P*dt
end where

end function ETD_uncoupled


function grad3D1_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 1-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE - EXCEPT FOR THE GRID.
!-------STRUCTURE.  
!-------
!-------FOR THIS AND FOLLOWING CURV DERIVATIVE PROCEDURES THE UPPER
!-------AND LOWER BOUNDS FOR DIFFERENTIATION CANNOT BE DETERMINED
!-------WITHOUT ADDITIONAL ARGUMENTS SINCE STRUCTURE X CANNOT BE
!-------BE TRIMMED IN THE SAME WAY AS A NORMAL ARRAY.  ESSENTIALLY
!-------THE INPUT ARRAY "F" AND GRID NEED TO BE INDEXED DIFFERENTLY, LEADING
!-------TO TWO SEPARATE SETS ON INDICES THROUGHOUT THIS AND SIMILAR ROUTINES.
!-------
!-------ONE ISSUE ADDRESSED HERE IS THAT THE METRIC FACTORS NEED TO KNOW
!-------WHAT PART OF THE GRID THAT THEY ARE BEING USED OVER...  IE
!-------IF ANY OF THESE IS FULL GRID THEN THE HALL VERSIONS SHOUDL
!-------SHOULD BE USED FOR THE METRIC FACTORS...
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3    !upper and lower bounds on the mesh indices for this derivative

integer :: ix2,ix3,lx1,lx2,lx3

!    real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: h1,h2,h3   !local references to the metric factors to be used in the derivative
!    real(wp), dimension(1:size(f,1)) :: dx1    !local reference to the backward difference
real(wp), dimension(:,:,:), pointer :: h1   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D1_curv_3


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad3D1 gradient function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not gauranteed to be contiguous in memory so I'm not sure this is the way to go
if (lx3<=x%lx3+4) then
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information'
else
  error stop '!!!  Array size is larger full mesh.'
end if
dx1=>x%dx1(lbnd1:ubnd1)


!NOW EXECUTE THE FINITE DIFFERENCES - NOTE THAT LOOP INDICES ARE MEANT TO INDEX ARRAY BEING DIFFERENCED AND NOT THE MESH STRUCTURE, WHICH USES INPUT BOUNDS.  TO KEEP THE CODE CLEAN I'VE ALIASED THE GRID VARS SO THAT THEY MAY BE ACCESSED BY LOOP INDEX.
do ix3=1,lx3
  do ix2=1,lx2
    grad3D1_curv_3(1,ix2,ix3)=(f(2,ix2,ix3)-f(1,ix2,ix3))/dx1(1)/h1(1,ix2,ix3)    !fwd diff. at beginning, note that h1 is cell-centered
    grad3D1_curv_3(2:lx1-1,ix2,ix3)=(f(3:lx1,ix2,ix3)-f(1:lx1-2,ix2,ix3)) &
                             /(dx1(3:lx1)+dx1(2:lx1-1))/h1(2:lx1-1,ix2,ix3)      !centered diff. in the middleq
    grad3D1_curv_3(lx1,ix2,ix3)=(f(lx1,ix2,ix3)-f(lx1-1,ix2,ix3))/dx1(lx1)/h1(lx1,ix2,ix3)    !backward diff. at end
  end do
end do
end function grad3D1_curv_3


function grad3D1_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 1-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE - EXCEPT FOR THE GRID.
!-------STRUCTURE.  
!-------
!-------FOR THIS AND FOLLOWING CURV DERIVATIVE PROCEDURES THE UPPER
!-------AND LOWER BOUNDS FOR DIFFERENTIATION CANNOT BE DETERMINED
!-------WITHOUT ADDITIONAL ARGUMENTS SINCE STRUCTURE X CANNOT BE
!-------BE TRIMMED IN THE SAME WAY AS A NORMAL ARRAY.  ESSENTIALLY
!-------THE INPUT ARRAY "F" AND GRID NEED TO BE INDEXED DIFFERENTLY, LEADING
!-------TO TWO SEPARATE SETS ON INDICES THROUGHOUT THIS AND SIMILAR ROUTINES.
!-------
!-------ONE ISSUE ADDRESSED HERE IS THAT THE METRIC FACTORS NEED TO KNOW
!-------WHAT PART OF THE GRID THAT THEY ARE BEING USED OVER...  IE
!-------IF ANY OF THESE IS FULL GRID THEN THE HALL VERSIONS SHOUDL
!-------SHOULD BE USED FOR THE METRIC FACTORS...
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3    !upper and lower bounds on the mesh indices for this derivative

integer :: ix2,ix3,lx1,lx2,lx3

!    real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: h1,h2,h3   !local references to the metric factors to be used in the derivative
!    real(wp), dimension(1:size(f,1)) :: dx1    !local reference to the backward difference
real(wp), dimension(:,:,:), pointer :: h1   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D1_curv_23


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad3D1 gradient function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not gauranteed to be contiguous in memory so I'm not sure this is the way to go
if (lx3<=x%lx3+4) then
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information'
else
  error stop '!!!  Array size is larger full mesh.'
end if
dx1=>x%dx1(lbnd1:ubnd1)


!NOW EXECUTE THE FINITE DIFFERENCES - NOTE THAT LOOP INDICES ARE MEANT TO INDEX ARRAY BEING DIFFERENCED AND NOT THE MESH STRUCTURE, WHICH USES INPUT BOUNDS.  TO KEEP THE CODE CLEAN I'VE ALIASED THE GRID VARS SO THAT THEY MAY BE ACCESSED BY LOOP INDEX.
do ix3=1,lx3
  do ix2=1,lx2
    grad3D1_curv_23(1,ix2,ix3)=(f(2,ix2,ix3)-f(1,ix2,ix3))/dx1(1)/h1(1,ix2,ix3)    !fwd diff. at beginning, note that h1 is cell-centered
    grad3D1_curv_23(2:lx1-1,ix2,ix3)=(f(3:lx1,ix2,ix3)-f(1:lx1-2,ix2,ix3)) &
                             /(dx1(3:lx1)+dx1(2:lx1-1))/h1(2:lx1-1,ix2,ix3)      !centered diff. in the middleq
    grad3D1_curv_23(lx1,ix2,ix3)=(f(lx1,ix2,ix3)-f(lx1-1,ix2,ix3))/dx1(lx1)/h1(lx1,ix2,ix3)    !backward diff. at end
  end do
end do
end function grad3D1_curv_23


function grad3D2_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h2   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx2    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D2_curv_3


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)

if (lx2>1) then    !if we have a singleton dimension then we are doing a 2D run and the derivatives in this direction are zero

  !ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
  if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
    error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
  end if


  !CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
  !INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
  if (lx3<=x%lx3+4) then
    h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  else if (lx3<=x%lx3all+4) then
    h2=>x%h2all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    print *,   '! Accessing root-only grid information'
  else
    error stop '!!!  Array size is larger full mesh.'
  end if
  dx2=>x%dx2(lbnd2:ubnd2)


  !DIFFERENCING
  do ix3=1,lx3
    do ix1=1,lx1
      grad3D2_curv_3(ix1,1,ix3)=(f(ix1,2,ix3)-f(ix1,1,ix3))/dx2(2)/h2(ix1,1,ix3)
      grad3D2_curv_3(ix1,2:lx2-1,ix3)=(f(ix1,3:lx2,ix3)-f(ix1,1:lx2-2,ix3)) &
                               /(dx2(3:lx2)+dx2(2:lx2-1))/h2(ix1,2:lx2-1,ix3)
      grad3D2_curv_3(ix1,lx2,ix3)=(f(ix1,lx2,ix3)-f(ix1,lx2-1,ix3))/dx2(lx2)/h2(ix1,lx2,lx3)
    end do
  end do
else
  grad3D2_curv_3=0._wp
end if
end function grad3D2_curv_3


function grad3D2_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h2   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx2    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D2_curv_23


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)

if (lx2>1) then    !if we have a singleton dimension then we are doing a 2D run and the derivatives in this direction are zero

  !ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
  if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
    error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
  end if


  !CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
  !INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
  if (lx3<=x%lx3+4) then     !this is a derivative over a slab region (subgrid)
    h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    dx2=>x%dx2(lbnd2:ubnd2)
  else if (lx3<=x%lx3all+4) then    !if a larger dimension was specified for x3 then assume that we are differentiating over x2all and x3all
    h2=>x%h2all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    dx2=>x%dx2all(lbnd2:ubnd2)
    print *,   '! Accessing root-only grid information, while taking derivative in 2-direction'
  else
    error stop '!!!  Array size is larger full mesh.'
  end if


  !DIFFERENCING
  do ix3=1,lx3
    do ix1=1,lx1
      grad3D2_curv_23(ix1,1,ix3)=(f(ix1,2,ix3)-f(ix1,1,ix3))/dx2(2)/h2(ix1,1,ix3)
      grad3D2_curv_23(ix1,2:lx2-1,ix3)=(f(ix1,3:lx2,ix3)-f(ix1,1:lx2-2,ix3)) &
                               /(dx2(3:lx2)+dx2(2:lx2-1))/h2(ix1,2:lx2-1,ix3)
      grad3D2_curv_23(ix1,lx2,ix3)=(f(ix1,lx2,ix3)-f(ix1,lx2-1,ix3))/dx2(lx2)/h2(ix1,lx2,lx3)
    end do
  end do
else
  grad3D2_curv_23=0._wp
end if
end function grad3D2_curv_23


function grad3D3_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x 
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix2,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D3_curv_3


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not gauranteed to be contiguous in memory so I'm not sure this is the way to go
if (lx3<=x%lx3+4) then
  h3=>x%h3(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3(lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then
  h3=>x%h3all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3all(lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


!FINITE DIFFERENCING
do ix2=1,lx2
  do ix1=1,lx1
    grad3D3_curv_3(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,1))/dx3(2)/h3(ix1,ix2,1)
    grad3D3_curv_3(ix1,ix2,2:lx3-1)=(f(ix1,ix2,3:lx3)-f(ix1,ix2,1:lx3-2)) &
                             /(dx3(3:lx3)+dx3(2:lx3-1))/h3(ix1,ix2,2:lx3-1)
    grad3D3_curv_3(ix1,ix2,lx3)=(f(ix1,ix2,lx3)-f(ix1,ix2,lx3-1))/dx3(lx3)/h3(ix1,ix2,lx3)
  end do
end do

end function grad3D3_curv_3


function grad3D3_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x 
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix2,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D3_curv_23


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


if (lx3>1) then    !only differentiate if we have non-singleton dimension, otherwise set to zero
  !ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
  if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
    error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
  end if
  
  
  !CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
  !INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
  !Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
  !Unfortunately a pointer is not gauranteed to be contiguous in memory so I'm not sure this is the way to go
  if (lx3<=x%lx3+4) then
    h3=>x%h3(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    dx3=>x%dx3(lbnd3:ubnd3)
  else if (lx3<=x%lx3all+4) then
    h3=>x%h3all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    dx3=>x%dx3all(lbnd3:ubnd3)
    print *,   '! Accessing root-only grid information while differentiating in x3'
  else
    error stop '!!!  Array size is larger than full mesh.'
  end if
  
  
  !FINITE DIFFERENCING
  do ix2=1,lx2
    do ix1=1,lx1
      grad3D3_curv_23(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,1))/dx3(2)/h3(ix1,ix2,1)
      grad3D3_curv_23(ix1,ix2,2:lx3-1)=(f(ix1,ix2,3:lx3)-f(ix1,ix2,1:lx3-2)) &
                               /(dx3(3:lx3)+dx3(2:lx3-1))/h3(ix1,ix2,2:lx3-1)
      grad3D3_curv_23(ix1,ix2,lx3)=(f(ix1,ix2,lx3)-f(ix1,ix2,lx3-1))/dx3(lx3)/h3(ix1,ix2,lx3)
    end do
  end do
else
  grad3D3_curv_23=0._wp
end if

end function grad3D3_curv_23


function div3D_curv_3(f1,f2,f3,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D DIVERGENCE.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  DX(I) IS PRESUMED
!-------TO BE THE *BACKWARD* DIFFERENCE AT POINT I
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f1,f2,f3
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix2,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h1,h2,h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1    
real(wp), dimension(:), pointer :: dx2
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

real(wp), dimension(1:size(f1,1),1:size(f1,2),1:size(f1,3)) :: div3D_curv_3


lx1=size(f1,1)
lx2=size(f1,2)
lx3=size(f1,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not gauranteed to be contiguous in memory so there may be a performance hit in doing this...
dx1=>x%dx1(lbnd1:ubnd1)
dx2=>x%dx2(lbnd2:ubnd2)
if (lx3<=x%lx3+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3(lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then     !presumes root.  may cause a seg fault since workers don't have the full-grid metric factors (should they?)
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3all(lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information in divergence function div3D'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


!REASSIGNING POINTERS AS SUBARRAYS TO BE DIFFERENCED PRODUCED CLEAN CODE (BELOW), BUT NOW REQUIRES THAT WE OFFSET
!ARRAY INDICES FOR THE POINTER VARIABLES, SINCE, BY DEFAULT, THEY START AT 1...  ACTUALLY THIS SHOUDL NOT BE THE CASE
!SINCE EVERYTHING IS GETTING OFFSET (F UPON BEING PASSED INTO THIS FUNCTION AND THEN DIFFERENTIALS IN THE POINTER
!ASSIGNMENTS...


!FINITE DIFFERENCES
do ix3=1,lx3
  do ix2=1,lx2
    div3D_curv_3(1,ix2,ix3)=(h2(2,ix2,ix3)*h3(2,ix2,ix3)*f1(2,ix2,ix3)-h2(1,ix2,ix3)*h3(1,ix2,ix3)*f1(1,ix2,ix3))/dx1(2)
    div3D_curv_3(2:lx1-1,ix2,ix3)=(h2(3:lx1,ix2,ix3)*h3(3:lx1,ix2,ix3)*f1(3:lx1,ix2,ix3)- &
                h2(1:lx1-2,ix2,ix3)*h3(1:lx1-2,ix2,ix3)*f1(1:lx1-2,ix2,ix3)) / (dx1(3:lx1)+dx1(2:lx1-1))
    div3D_curv_3(lx1,ix2,ix3)=(h2(lx1,ix2,ix3)*h3(lx1,ix2,ix3)*f1(lx1,ix2,ix3)- &
                             h2(lx1-1,ix2,ix3)*h3(lx1-1,ix2,ix3)*f1(lx1-1,ix2,ix3))/dx1(lx1)
  end do
end do

if (lx2>1) then   !only if the x2-direction is not null
  do ix3=1,lx3
    do ix1=1,lx1
      div3D_curv_3(ix1,1,ix3)=div3D_curv_3(ix1,1,ix3)+ &
                    (h1(ix1,2,ix3)*h3(ix1,2,ix3)*f2(ix1,2,ix3)-h1(ix1,1,ix3)*h3(ix1,1,ix3)*f2(ix1,1,ix3))/dx2(2)
      div3D_curv_3(ix1,2:lx2-1,ix3)=div3D_curv_3(ix1,2:lx2-1,ix3)+ & 
                    (h1(ix1,3:lx2,ix3)*h3(ix1,3:lx2,ix3)*f2(ix1,3:lx2,ix3)- &
                     h1(ix1,1:lx2-2,ix3)*h3(ix1,1:lx2-2,ix3)*f2(ix1,1:lx2-2,ix3)) &
                             /(dx2(3:lx2)+dx2(2:lx2-1))
      div3D_curv_3(ix1,lx2,ix3)=div3D_curv_3(ix1,lx2,ix3)+ &
                    (h1(ix1,lx2,ix3)*h3(ix1,lx2,ix3)*f2(ix1,lx2,ix3)- &
                     h1(ix1,lx2-1,ix3)*h3(ix1,lx2-1,ix3)*f2(ix1,lx2-1,ix3))/dx2(lx2)
    end do
  end do
end if

do ix2=1,lx2
  do ix1=1,lx1
    div3D_curv_3(ix1,ix2,1)=div3D_curv_3(ix1,ix2,1)+ &
               (h1(ix1,ix2,2)*h2(ix1,ix2,2)*f3(ix1,ix2,2)-h1(ix1,ix2,1)*h2(ix1,ix2,1)*f3(ix1,ix2,1))/dx3(2)
    div3D_curv_3(ix1,ix2,2:lx3-1)=div3D_curv_3(ix1,ix2,2:lx3-1)+ &
               (h1(ix1,ix2,3:lx3)*h2(ix1,ix2,3:lx3)*f3(ix1,ix2,3:lx3)- &
                h1(ix1,ix2,1:lx3-2)*h2(ix1,ix2,1:lx3-2)*f3(ix1,ix2,1:lx3-2))&
                          /(dx3(3:lx3)+dx3(2:lx3-1))
    div3D_curv_3(ix1,ix2,lx3)=div3D_curv_3(ix1,ix2,lx3)+ &
               (h1(ix1,ix2,lx3)*h2(ix1,ix2,lx3)*f3(ix1,ix2,lx3)- &
                h1(ix1,ix2,lx3-1)*h2(ix1,ix2,lx3-1)*f3(ix1,ix2,lx3-1))/dx3(lx3)
  end do
end do

div3D_curv_3=div3D_curv_3/(h1*h2*h3)

end function div3D_curv_3


function div3D_curv_23(f1,f2,f3,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)

!------------------------------------------------------------
!-------COMPUTE A 3D DIVERGENCE.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  DX(I) IS PRESUMED
!-------TO BE THE *BACKWARD* DIFFERENCE AT POINT I
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f1,f2,f3
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3

integer :: ix1,ix2,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h1,h2,h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1    
real(wp), dimension(:), pointer :: dx2
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

real(wp), dimension(1:size(f1,1),1:size(f1,2),1:size(f1,3)) :: div3D_curv_23


lx1=size(f1,1)
lx2=size(f1,2)
lx3=size(f1,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not gauranteed to be contiguous in memory so there may be a performance hit in doing this...
dx1=>x%dx1(lbnd1:ubnd1)
if (lx3<=x%lx3+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3(lbnd3:ubnd3)
  dx2=>x%dx2(lbnd2:ubnd2)
else if (lx3<=x%lx3all+4) then     !presumes root or some process that has access to ALL full grid variables (normally only root).
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3all(lbnd3:ubnd3)
  dx2=>x%dx2all(lbnd2:ubnd2)
  print *,   '! Accessing root-only grid information in divergence function div3D'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


!REASSIGNING POINTERS AS SUBARRAYS TO BE DIFFERENCED PRODUCED CLEAN CODE (BELOW), BUT NOW REQUIRES THAT WE OFFSET
!ARRAY INDICES FOR THE POINTER VARIABLES, SINCE, BY DEFAULT, THEY START AT 1...  ACTUALLY THIS SHOUDL NOT BE THE CASE
!SINCE EVERYTHING IS GETTING OFFSET (F UPON BEING PASSED INTO THIS FUNCTION AND THEN DIFFERENTIALS IN THE POINTER
!ASSIGNMENTS...


!FINITE DIFFERENCES
do ix3=1,lx3
  do ix2=1,lx2
    div3D_curv_23(1,ix2,ix3)=(h2(2,ix2,ix3)*h3(2,ix2,ix3)*f1(2,ix2,ix3)-h2(1,ix2,ix3)*h3(1,ix2,ix3)*f1(1,ix2,ix3))/dx1(2)
    div3D_curv_23(2:lx1-1,ix2,ix3)=(h2(3:lx1,ix2,ix3)*h3(3:lx1,ix2,ix3)*f1(3:lx1,ix2,ix3)- &
                h2(1:lx1-2,ix2,ix3)*h3(1:lx1-2,ix2,ix3)*f1(1:lx1-2,ix2,ix3)) / (dx1(3:lx1)+dx1(2:lx1-1))
    div3D_curv_23(lx1,ix2,ix3)=(h2(lx1,ix2,ix3)*h3(lx1,ix2,ix3)*f1(lx1,ix2,ix3)- &
                             h2(lx1-1,ix2,ix3)*h3(lx1-1,ix2,ix3)*f1(lx1-1,ix2,ix3))/dx1(lx1)
  end do
end do

if (lx2>1) then   !only if the x2-direction is not null
  do ix3=1,lx3
    do ix1=1,lx1
      div3D_curv_23(ix1,1,ix3)=div3D_curv_23(ix1,1,ix3)+ &
                    (h1(ix1,2,ix3)*h3(ix1,2,ix3)*f2(ix1,2,ix3)-h1(ix1,1,ix3)*h3(ix1,1,ix3)*f2(ix1,1,ix3))/dx2(2)
      div3D_curv_23(ix1,2:lx2-1,ix3)=div3D_curv_23(ix1,2:lx2-1,ix3)+ & 
                    (h1(ix1,3:lx2,ix3)*h3(ix1,3:lx2,ix3)*f2(ix1,3:lx2,ix3)- &
                     h1(ix1,1:lx2-2,ix3)*h3(ix1,1:lx2-2,ix3)*f2(ix1,1:lx2-2,ix3)) &
                             /(dx2(3:lx2)+dx2(2:lx2-1))
      div3D_curv_23(ix1,lx2,ix3)=div3D_curv_23(ix1,lx2,ix3)+ &
                    (h1(ix1,lx2,ix3)*h3(ix1,lx2,ix3)*f2(ix1,lx2,ix3)- &
                     h1(ix1,lx2-1,ix3)*h3(ix1,lx2-1,ix3)*f2(ix1,lx2-1,ix3))/dx2(lx2)
    end do
  end do
end if

if (lx3>1) then    !only if non-singleton 3rd dimension
  do ix2=1,lx2
    do ix1=1,lx1
      div3D_curv_23(ix1,ix2,1)=div3D_curv_23(ix1,ix2,1)+ &
                 (h1(ix1,ix2,2)*h2(ix1,ix2,2)*f3(ix1,ix2,2)-h1(ix1,ix2,1)*h2(ix1,ix2,1)*f3(ix1,ix2,1))/dx3(2)
      div3D_curv_23(ix1,ix2,2:lx3-1)=div3D_curv_23(ix1,ix2,2:lx3-1)+ &
                 (h1(ix1,ix2,3:lx3)*h2(ix1,ix2,3:lx3)*f3(ix1,ix2,3:lx3)- &
                  h1(ix1,ix2,1:lx3-2)*h2(ix1,ix2,1:lx3-2)*f3(ix1,ix2,1:lx3-2))&
                            /(dx3(3:lx3)+dx3(2:lx3-1))
      div3D_curv_23(ix1,ix2,lx3)=div3D_curv_23(ix1,ix2,lx3)+ &
                 (h1(ix1,ix2,lx3)*h2(ix1,ix2,lx3)*f3(ix1,ix2,lx3)- &
                  h1(ix1,ix2,lx3-1)*h2(ix1,ix2,lx3-1)*f3(ix1,ix2,lx3-1))/dx3(lx3)
    end do
  end do
end if

div3D_curv_23=div3D_curv_23/(h1*h2*h3)

end function div3D_curv_23



!THE REMAINDER OF THE FUNCTIONS IN THIS MODULE DO *NOT* TAKE METRIC FACTORS INTO ACCOUNT.
!HENCE THEY ARE PURE DERIVATIVES AND INTEGRALS...  These probably need to be renamed and
!also some error checking should probably be done...

pure function grad2D1_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 1-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix2,lx1,lx2

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D1_curv

lx1=size(f,1)
lx2=size(f,2)

do ix2=1,lx2
  grad2D1_curv(1,ix2)=(f(2,ix2)-f(1,ix2))/x%dx1(lbnd+1)
  grad2D1_curv(2:lx1-1,ix2)=(f(3:lx1,ix2)-f(1:lx1-2,ix2)) &
                           /(x%dx1(lbnd+2:ubnd)+x%dx1(lbnd+1:ubnd-1))
  grad2D1_curv(lx1,ix2)=(f(lx1,ix2)-f(lx1-1,ix2))/x%dx1(ubnd)
end do
end function grad2D1_curv


pure function grad2D1_curv_alt_3(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 1-DIMENSION USING THE
!-------VARIABLE X2 AS THE DIFFERENTIAL. IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix2,lx1,lx2

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D1_curv_alt_3

lx1=size(f,1)
lx2=size(f,2)

do ix2=1,lx2
  grad2D1_curv_alt_3(1,ix2)=(f(2,ix2)-f(1,ix2))/x%dx2(lbnd+1)
  grad2D1_curv_alt_3(2:lx1-1,ix2)=(f(3:lx1,ix2)-f(1:lx1-2,ix2)) &
                           /(x%dx2(lbnd+2:ubnd)+x%dx2(lbnd+1:ubnd-1))
  grad2D1_curv_alt_3(lx1,ix2)=(f(lx1,ix2)-f(lx1-1,ix2))/x%dx2(ubnd)
end do
end function grad2D1_curv_alt_3


function grad2D1_curv_alt_23(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 1-DIMENSION USING THE
!-------VARIABLE X2 AS THE DIFFERENTIAL. IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

!may need to explicitly check that we aren't differentiating over lx2=1???

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix3,lx2,lx3
real(wp), dimension(:), pointer :: dx2

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D1_curv_alt_23

lx2=size(f,1)
lx3=size(f,2)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx2 /= ubnd-lbnd+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad2D1_curv_alt_23.'   !just bail on it and let the user figure it out
end if


if (lx2<=x%lx2+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  dx2=>x%dx2(lbnd:ubnd)
else if (lx2<=x%lx2all+4) then     !presumes root or some process that has access to ALL full grid variables (normally only root).
  dx2=>x%dx2all(lbnd:ubnd)
  print *,   '! Accessing root-only grid information in divergence function grad2D1_curv_alt_23'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


!What follows is slightly confusing because the 2-dimension is partly indexed by lx1
do ix3=1,lx3
  grad2D1_curv_alt_23(1,ix3)=(f(2,ix3)-f(1,ix3))/dx2(2)
  grad2D1_curv_alt_23(2:lx2-1,ix3)=(f(3:lx2,ix3)-f(1:lx2-2,ix3)) &
                           /(dx2(3:lx2)+dx2(2:lx2-1))
  grad2D1_curv_alt_23(lx2,ix3)=(f(lx2,ix3)-f(lx2-1,ix3))/dx2(lx2)
end do
end function grad2D1_curv_alt_23


pure function grad2D2_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,lx1,lx2

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D2_curv

lx1=size(f,1)
lx2=size(f,2)

if (lx2>1) then         !are we even simulating x2-direction?
  do ix1=1,lx1
    grad2D2_curv(ix1,1)=(f(ix1,2)-f(ix1,1))/x%dx2(lbnd+1)
    grad2D2_curv(ix1,2:lx2-1)=(f(ix1,3:lx2)-f(ix1,1:lx2-2)) &
                             /(x%dx2(lbnd+2:ubnd)+x%dx2(lbnd+1:ubnd-1))
    grad2D2_curv(ix1,lx2)=(f(ix1,lx2)-f(ix1,lx2-1))/x%dx2(ubnd)
  end do
else
  grad2D2_curv=0._wp
end if
end function grad2D2_curv


function grad2D2_curv_23(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,lx1,lx2
real(wp), dimension(:), pointer :: dx2

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D2_curv_23

lx1=size(f,1)
lx2=size(f,2)

!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx2 /= ubnd-lbnd+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad2D2_curv_23.'   !just bail on it and let the user figure it out
end if


if (lx2<=x%lx2+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  dx2=>x%dx2(lbnd:ubnd)
else if (lx2<=x%lx2all+4) then     !presumes root or some process that has access to ALL full grid variables (normally only root).
  dx2=>x%dx2all(lbnd:ubnd)
  print *,   '! Accessing root-only grid information in divergence function grad2D2_curv_23'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


if (lx2>1) then         !are we even simulating x2-direction?
  do ix1=1,lx1
    grad2D2_curv_23(ix1,1)=(f(ix1,2)-f(ix1,1))/dx2(2)
    grad2D2_curv_23(ix1,2:lx2-1)=(f(ix1,3:lx2)-f(ix1,1:lx2-2)) &
                             /(dx2(3:lx2)+dx2(2:lx2-1))
    grad2D2_curv_23(ix1,lx2)=(f(ix1,lx2)-f(ix1,lx2-1))/dx2(lx2)
  end do
else
  grad2D2_curv_23=0._wp
end if
end function grad2D2_curv_23


pure function grad2D3_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  A SEPARATE ROUTINE
!-------IS NEEDED FOR THE CURVILINEAR CASE, BECAUSE OTHERWISE
!-------WE HAVE NO WAY TO KNOW WHAT DIMENSION IN THE STRUCTURE
!-------SHOULD BE USED FOR THE DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,lx1,lx3

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D3_curv

lx1=size(f,1)    !this is really the 2-dimension for a flattened array
lx3=size(f,2)


if (lx3 == x%lx3all) then    !derivative over full grid (could maybe use some error checking with)
  do ix1=1,lx1
    grad2D3_curv(ix1,1)=(f(ix1,2)-f(ix1,1))/x%dx3all(lbnd+1)
    grad2D3_curv(ix1,2:lx3-1)=(f(ix1,3:lx3)-f(ix1,1:lx3-2)) &
                             /(x%dx3all(lbnd+2:ubnd)+x%dx3all(lbnd+1:ubnd-1))
    grad2D3_curv(ix1,lx3)=(f(ix1,lx3)-f(ix1,lx3-1))/x%dx3all(ubnd)
  end do
else                         !derivative over part of the grid
  do ix1=1,lx1
    grad2D3_curv(ix1,1)=(f(ix1,2)-f(ix1,1))/x%dx3(lbnd+1)
    grad2D3_curv(ix1,2:lx3-1)=(f(ix1,3:lx3)-f(ix1,1:lx3-2)) &
                             /(x%dx3(lbnd+2:ubnd)+x%dx3(lbnd+1:ubnd-1))
    grad2D3_curv(ix1,lx3)=(f(ix1,lx3)-f(ix1,lx3-1))/x%dx3(ubnd)
  end do
end if
end function grad2D3_curv


function grad2D3_curv_23(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  A SEPARATE ROUTINE
!-------IS NEEDED FOR THE CURVILINEAR CASE, BECAUSE OTHERWISE
!-------WE HAVE NO WAY TO KNOW WHAT DIMENSION IN THE STRUCTURE
!-------SHOULD BE USED FOR THE DERIVATIVE.
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix2,lx2,lx3
real(wp), dimension(:), pointer :: dx3

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D3_curv_23

lx2=size(f,1)    !this is really the 2-dimension for a flattened array
lx3=size(f,2)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx3 /= ubnd-lbnd+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad2D3_curv_alt_23.'   !just bail on it and let the user figure it out
end if


if (lx3<=x%lx3+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  dx3=>x%dx3(lbnd:ubnd)
else if (lx3<=x%lx3all+4) then     !presumes root or some process that has access to ALL full grid variables (normally only root).
  dx3=>x%dx3all(lbnd:ubnd)
  print *,   '! Accessing root-only grid information in divergence function grad2D3_curv_alt_23'
else
  error stop '!!!  Array size is larger than full mesh.'
end if


do ix2=1,lx2
  grad2D3_curv_23(ix2,1)=(f(ix2,2)-f(ix2,1))/dx3(2)
  grad2D3_curv_23(ix2,2:lx3-1)=(f(ix2,3:lx3)-f(ix2,1:lx3-2)) &
                           /(dx3(3:lx3)+dx3(2:lx3-1))
  grad2D3_curv_23(ix2,lx3)=(f(ix2,lx3)-f(ix2,lx3-1))/dx3(lx3)
end do

end function grad2D3_curv_23


pure function integral3D1_curv(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:512:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral3D1_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 3D ARRAY ALONG THE 1-DIM.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------NOTE, AS BEFORE THAT THE GRID MAY BE INDEXED DIFFERENTLY
!-------THAT THE FUNCTION BEING INTEGRATED
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd    !no upper bound used, but could be for error checking

integer :: ix1,lx1
real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: integral3D1_curv

lx1=size(f,1)

integral3D1_curv(1,:,:)=0._wp
do ix1=2,lx1
  integral3D1_curv(ix1,:,:)=integral3D1_curv(ix1-1,:,:)+0.5d0*(f(ix1,:,:)+f(ix1-1,:,:))*x%dx1(lbnd+ix1-1)
end do
end function integral3D1_curv


pure function integral3D1_curv_alt(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:539:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral3D1_curv_alt(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 3D ARRAY ALONG THE 1-DIM.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------NOTE, AS BEFORE THAT THE GRID MAY BE INDEXED DIFFERENTLY
!-------THAN THE FUNCTION BEING INTEGRATED.  ALSO THIS ALT VERSION
!-------INTEGRATES FROM SET OPINT TO THE MAX VALUE ON THE GRID, 
!-------RATHER THAN FROM THE MIN VALUE TO A FIXED POINT.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd    !no upper bound used, but could be for error checking

integer :: ix1,lx1
real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: integral3D1_curv_alt

lx1=size(f,1)

integral3D1_curv_alt(lx1,:,:)=0._wp
do ix1=lx1-1,1,-1    !start from the logical top and sum  downward (keep dx positive since the intent is to integrate from fixed point to top)
  integral3D1_curv_alt(ix1,:,:)=integral3D1_curv_alt(ix1+1,:,:)+0.5d0*(f(ix1,:,:)+f(ix1+1,:,:))*x%dx1(lbnd+ix1-1+1)   !+1 since we are starting from top
end do
end function integral3D1_curv_alt


pure function integral2D1_curv(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:568:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral2D1_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 2D ARRAY ALONG THE 1-DIM.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,lx1
real(wp), dimension(1:size(f,1),1:size(f,2)) :: integral2D1_curv

lx1=size(f,1)

integral2D1_curv(1,:)=0._wp
do ix1=2,lx1
  integral2D1_curv(ix1,:)=integral2D1_curv(ix1-1,:)+0.5d0*(f(ix1,:)+f(ix1-1,:))*x%dx1(lbnd+ix1-1)
end do
end function integral2D1_curv


pure function integral2D1_curv_alt(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:592:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral2D1_curv_alt(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 2D ARRAY ALONG THE 1-DIM.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.
!-------
!-------THIS VERSION INTEGRATES WRT X2
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix2,lx2
real(wp), dimension(1:size(f,1),1:size(f,2)) :: integral2D1_curv_alt

lx2=size(f,1)

integral2D1_curv_alt(1,:)=0._wp
do ix2=2,lx2
  integral2D1_curv_alt(ix2,:)=integral2D1_curv_alt(ix2-1,:)+0.5d0*(f(ix2,:)+f(ix2-1,:))*x%dx2(lbnd+ix2-1)
end do
end function integral2D1_curv_alt


pure function integral2D2_curv(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:618:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral2D2_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 2D ARRAY ALONG THE 2-DIM.  
!-------WITH RESPECT TO X2.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------IT IS ASSUMED THAT THE INTEGRATION IS ALWAYS OVER THE
!-------ENTIRE GRID (X2ALL).
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix2,lx2
real(wp), dimension(1:size(f,1),1:size(f,2)) :: integral2D2_curv

lx2=size(f,2)

integral2D2_curv(:,1)=0._wp
do ix2=2,lx2
  integral2D2_curv(:,ix2)=integral2D2_curv(:,ix2-1)+0.5d0*(f(:,ix2)+f(:,ix2-1))*x%dx2all(lbnd+ix2-1)
end do
end function integral2D2_curv


pure function integral2D2_curv_alt(f,x,lbnd,ubnd)

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:642:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral2D2_curv_alt(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 2D ARRAY ALONG THE 2-DIM, 
!-------WITH RESPECT TO X3.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.
!-------
!-------IT IS ASSUMED THAT THE INTEGRATION IS ALWAYS OVER THE
!-------ENTIRE GRID (X3ALL).
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix3,lx3
real(wp), dimension(1:size(f,1),1:size(f,2)) :: integral2D2_curv_alt

lx3=size(f,2)

integral2D2_curv_alt(:,1)=0._wp
do ix3=2,lx3
  integral2D2_curv_alt(:,ix3)=integral2D2_curv_alt(:,ix3-1)+0.5d0*(f(:,ix3)+f(:,ix3-1))*x%dx3all(lbnd+ix3-1)
end do
end function integral2D2_curv_alt


pure function grad3D3_curv_periodic(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  THIS ROUTINE ASSUMES
!-------A PERIODIC BOOUNDARY IN THE X3 DIMENSION
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!-------
!-------THIS FUNCTION IS ONLY VALID WITH CARTESIAN GRIDS.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,ix2,lx1,lx2,lx3

real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: grad3D3_curv_periodic


lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)

if(lx3 == x%lx3all) then    !differentiating over full grid
  do ix2=1,lx2
    do ix1=1,lx1
!          grad3D3_curv_periodic(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,1))/x%dx3all(lbnd+1)
      grad3D3_curv_periodic(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,lx3))/(x%dx3all(lbnd+1)+x%dx3all(lbnd))    !this assumes that the backward difference for hte first cell has been set the the same as the forward difference for the final cell, i.e. x%dx3all(lbnd)==x%dx3all(ubnd+1).  In general when doing periodic grids it is probably best to hard code all of the differences outside the domain to be equal (or to use uniform meshes)
      grad3D3_curv_periodic(ix1,ix2,2:lx3-1)=(f(ix1,ix2,3:lx3)-f(ix1,ix2,1:lx3-2)) &
                               /(x%dx3all(lbnd+2:ubnd)+x%dx3all(lbnd+1:ubnd-1))
      grad3D3_curv_periodic(ix1,ix2,lx3)=(f(ix1,ix2,1)-f(ix1,ix2,lx3-1))/(x%dx3all(ubnd)+x%dx3all(ubnd+1))
    end do
  end do
else
  do ix2=1,lx2
    do ix1=1,lx1
!          grad3D3_curv_periodic(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,1))/x%dx3(lbnd+1)
      grad3D3_curv_periodic(ix1,ix2,1)=(f(ix1,ix2,2)-f(ix1,ix2,lx3))/(x%dx3(lbnd+1)+x%dx3(lbnd))
      grad3D3_curv_periodic(ix1,ix2,2:lx3-1)=(f(ix1,ix2,3:lx3)-f(ix1,ix2,1:lx3-2)) &
                               /(x%dx3(lbnd+2:ubnd)+x%dx3(lbnd+1:ubnd-1))
!          grad3D3_curv_periodic(ix1,ix2,lx3)=(f(ix1,ix2,lx3)-f(ix1,ix2,lx3-1))/x%dx3(ubnd)
      grad3D3_curv_periodic(ix1,ix2,lx3)=(f(ix1,ix2,1)-f(ix1,ix2,lx3-1))/(x%dx3(ubnd)+x%dx3(ubnd+1))
    end do
  end do
end if

end function grad3D3_curv_periodic


pure function grad2D3_curv_periodic(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE A 2D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT 
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE.  A SEPARATE ROUTINE
!-------IS NEEDED FOR THE CURVILINEAR CASE, BECAUSE OTHERWISE
!-------WE HAVE NO WAY TO KNOW WHAT DIMENSION IN THE STRUCTURE
!-------SHOULD BE USED FOR THE DERIVATIVE.  NOTE THAT THIS USES
!-------A PERIODIC DIFFERENCE IF THE DERIVATIVE IS OVER THE ENTIRE
!-------GRID AND APERIODIC OTHERWISE.
!-------
!-------THIS FUNCTION IS ONLY VALID WITH CARTESIAN GRIDS.
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: f
type(curvmesh), intent(in) :: x
integer, intent(in) :: lbnd,ubnd

integer :: ix1,lx1,lx3

real(wp), dimension(1:size(f,1),1:size(f,2)) :: grad2D3_curv_periodic

lx1=size(f,1)    !this is really the 2-dimension for a flattened array
lx3=size(f,2)

if (lx3 == x%lx3all) then    !derivative over full grid (could maybe use some error checking with)
  do ix1=1,lx1
!        grad2D3_curv_periodic(ix1,1)=(f(ix1,2)-f(ix1,1))/x%dx3all(lbnd+1)
    grad2D3_curv_periodic(ix1,1)=(f(ix1,2)-f(ix1,lx3))/(x%dx3all(lbnd+1)+x%dx3all(lbnd))    !this assumes that the backward difference for hte first cell has been set the the same as the forward difference for the final cell, i.e. x%dx3all(lbnd)==x%dx3all(ubnd+1).  In general when doing periodic grids it is probably best to hard code all of the differences outside the domain to be equal (or to use uniform meshes)
    grad2D3_curv_periodic(ix1,2:lx3-1)=(f(ix1,3:lx3)-f(ix1,1:lx3-2)) &
                             /(x%dx3all(lbnd+2:ubnd)+x%dx3all(lbnd+1:ubnd-1))
!        grad2D3_curv_periodic(ix1,lx3)=(f(ix1,lx3)-f(ix1,lx3-1))/x%dx3all(ubnd)
    grad2D3_curv_periodic(ix1,lx3)=(f(ix1,1)-f(ix1,lx3-1))/(x%dx3all(ubnd)+x%dx3all(ubnd+1))
  end do
else                         !derivative over part of the grid
  do ix1=1,lx1
    grad2D3_curv_periodic(ix1,1)=(f(ix1,2)-f(ix1,1))/x%dx3(lbnd+1)
    grad2D3_curv_periodic(ix1,2:lx3-1)=(f(ix1,3:lx3)-f(ix1,1:lx3-2)) &
                             /(x%dx3(lbnd+2:ubnd)+x%dx3(lbnd+1:ubnd-1))
    grad2D3_curv_periodic(ix1,lx3)=(f(ix1,lx3)-f(ix1,lx3-1))/x%dx3(ubnd)
  end do
end if
end function grad2D3_curv_periodic


pure function chapman_a(alt,nmax,alt0,H)

!------------------------------------------------------------
!-------GENERATE A CHAPMAN ALPHA PROFILE.
!------------------------------------------------------------


real(wp), dimension(:,:,:), intent(in) :: alt,H
real(wp), intent(in) :: nmax,alt0

real(wp), dimension(1:size(alt,1),1:size(alt,2),1:size(alt,3)) :: chapman_a

chapman_a=nmax*exp(0.5_wp*(1.0_wp-(alt-alt0)/H-exp(-(alt-alt0)/H)))

where (chapman_a<1.0_wp)
  chapman_a=1.0_wp
end where

end function chapman_a


end module calculus

