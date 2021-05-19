submodule (calculus) gradient

implicit none (type, external)

contains

module procedure grad3D1_curv_3
! grad3D1_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
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
!-------IF ANY OF THESE IS FULL GRID THEN THE HALL VERSIONS SHOULD
!-------SHOULD BE USED FOR THE METRIC FACTORS...
!------------------------------------------------------------

integer :: ix2,ix3,lx1,lx2,lx3

!    real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: h1,h2,h3
!! local references to the metric factors to be used in the derivative
!    real(wp), dimension(1:size(f,1)) :: dx1
!! local reference to the backward difference
real(wp), dimension(:,:,:), pointer :: h1
!! local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1
!! local reference to the backward difference

lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


!! ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad3D1 gradient function.'
  !! just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not guaranteed to be contiguous in memory so I'm not sure this is the way to go
if (lx3<=x%lx3+4) then
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information'
else
  error stop '!!!  Array size is larger full mesh.'
end if
dx1=>x%dx1(lbnd1:ubnd1)


!! NOW EXECUTE THE FINITE DIFFERENCES -
!! NOTE THAT LOOP INDICES ARE MEANT TO INDEX ARRAY BEING DIFFERENCED AND NOT THE MESH STRUCTURE,
!! WHICH USES INPUT BOUNDS.  TO KEEP THE CODE CLEAN I'VE ALIASED THE GRID VARS SO THAT THEY MAY BE ACCESSED BY LOOP INDEX.
do ix3=1,lx3
  do ix2=1,lx2
    grad3D1_curv_3(1,ix2,ix3) = (f(2,ix2,ix3)-f(1,ix2,ix3))/dx1(1)/h1(1,ix2,ix3)
    !! fwd diff. at beginning, note that h1 is cell-centered
    grad3D1_curv_3(2:lx1-1,ix2,ix3) = (f(3:lx1,ix2,ix3)-f(1:lx1-2,ix2,ix3)) &
                             /(dx1(3:lx1)+dx1(2:lx1-1))/h1(2:lx1-1,ix2,ix3)
    !! centered diff. in the middleq
    grad3D1_curv_3(lx1,ix2,ix3) = (f(lx1,ix2,ix3)-f(lx1-1,ix2,ix3))/dx1(lx1)/h1(lx1,ix2,ix3)
    !! backward diff. at end
  end do
end do
end procedure grad3D1_curv_3


module procedure grad3D1_curv_23
! grad3D1_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
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
!-------IF ANY OF THESE IS FULL GRID THEN THE HALL VERSIONS SHOULD
!-------SHOULD BE USED FOR THE METRIC FACTORS...
!------------------------------------------------------------

integer :: ix2,ix3,lx1,lx2,lx3

!    real(wp), dimension(1:size(f,1),1:size(f,2),1:size(f,3)) :: h1,h2,h3
!! local references to the metric factors to be used in the derivative
!    real(wp), dimension(1:size(f,1)) :: dx1
!! local reference to the backward difference
real(wp), dimension(:,:,:), pointer :: h1
!! local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1
!! local reference to the backward difference

lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in grad3D1 gradient function.'
  !! just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not guaranteed to be contiguous in memory so I'm not sure this is the way to go
if (lx3<=x%lx3+4) then
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
else if (lx3<=x%lx3all+4) then
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  print *,   '! Accessing root-only grid information'
else
  error stop '!!!  Array size is larger full mesh.'
end if
dx1=>x%dx1(lbnd1:ubnd1)


!! NOW EXECUTE THE FINITE DIFFERENCES
!! NOTE THAT LOOP INDICES ARE MEANT TO INDEX ARRAY BEING DIFFERENCED AND NOT THE MESH STRUCTURE,
!! WHICH USES INPUT BOUNDS.  TO KEEP THE CODE CLEAN I'VE ALIASED THE GRID VARS SO THAT THEY MAY BE ACCESSED BY LOOP INDEX.
do ix3=1,lx3
  do ix2=1,lx2
    grad3D1_curv_23(1,ix2,ix3) = (f(2,ix2,ix3)-f(1,ix2,ix3))/dx1(1)/h1(1,ix2,ix3)
    !! fwd diff. at beginning, note that h1 is cell-centered
    grad3D1_curv_23(2:lx1-1,ix2,ix3) = (f(3:lx1,ix2,ix3)-f(1:lx1-2,ix2,ix3)) &
                             /(dx1(3:lx1)+dx1(2:lx1-1))/h1(2:lx1-1,ix2,ix3)
    !! centered diff. in the middleq
    grad3D1_curv_23(lx1,ix2,ix3) = (f(lx1,ix2,ix3)-f(lx1-1,ix2,ix3))/dx1(lx1)/h1(lx1,ix2,ix3)
    !! backward diff. at end
  end do
end do
end procedure grad3D1_curv_23


module procedure grad3D2_curv_3
! grad3D2_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

integer :: ix1,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h2   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx2    !local reference to the backward difference

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
end procedure grad3D2_curv_3


module procedure grad3D2_curv_23
! grad3D2_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 2-DIMENSION.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------

integer :: ix1,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h2   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx2    !local reference to the backward difference

lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)

if (x%lx2>1) then    !if we have a singleton dimension then we are doing a 2D run and the derivatives in this direction are zero

  !ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
  if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
    error stop '!!!  Inconsistent array and mesh sizes in gradient function.'
    !! just bail on it and let the user figure it out
  end if


  !CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
  !INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
  if (lx3<=x%lx3+4) then     !this is a derivative over a slab region (subgrid)
    h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
    dx2=>x%dx2(lbnd2:ubnd2)
  else if (lx3<=x%lx3all+4) then
    !! if a larger dimension was specified for x3 then assume that we are differentiating over x2all and x3all
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
end procedure grad3D2_curv_23


module procedure grad3D3_curv_3
! grad3D3_curv_3(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

integer :: ix1,ix2,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

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
!Unfortunately a pointer is not guaranteed to be contiguous in memory so I'm not sure this is the way to go
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

end procedure grad3D3_curv_3


module procedure grad3D3_curv_23
! grad3D3_curv_23(f,x,lbnd1,ubnd1,lbnd2,ubnd2,lbnd3,ubnd3)
!------------------------------------------------------------
!-------COMPUTE A 3D GRADIENT ALONG THE 3-DIMENSION.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!-------
!-------AN EXTRA STEP IS NEEDED IN THIS ROUTINE SINCE WE HAVE
!-------TO DETERMINE WHETHER THIS IS A FULL-GRID OR SUBGRID
!-------DERIVATIVE.
!------------------------------------------------------------

integer :: ix1,ix2,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

lx1=size(f,1)
lx2=size(f,2)
lx3=size(f,3)


if (x%lx3>1) then    !only differentiate if we have non-singleton dimension, otherwise set to zero
  !ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
  if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
    error stop '!!!  Inconsistent array and mesh sizes in gradient function.'   !just bail on it and let the user figure it out
  end if


  !CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
  !INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
  !Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
  !Unfortunately a pointer is not guaranteed to be contiguous in memory so I'm not sure this is the way to go
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

end procedure grad3D3_curv_23

end submodule gradient
