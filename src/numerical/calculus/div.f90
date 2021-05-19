submodule (calculus) div

implicit none (type, external)

contains

module procedure div3D_curv_3

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

integer :: ix1,ix2,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h1,h2,h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1
real(wp), dimension(:), pointer :: dx2
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference


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
!Unfortunately a pointer is not guaranteed to be contiguous in memory so there may be a performance hit in doing this...
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
!ARRAY INDICES FOR THE POINTER VARIABLES, SINCE, BY DEFAULT, THEY START AT 1...  ACTUALLY THIS SHOULD NOT BE THE CASE
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

end procedure div3D_curv_3


module procedure div3D_curv_23

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

integer :: ix1,ix2,ix3,lx1,lx2,lx3

real(wp), dimension(:,:,:), pointer :: h1,h2,h3   !local references to the metric factors to be used in the derivative
real(wp), dimension(:), pointer :: dx1
real(wp), dimension(:), pointer :: dx2
real(wp), dimension(:), pointer :: dx3    !local reference to the backward difference

lx1=size(f1,1)
lx2=size(f1,2)
lx3=size(f1,3)


if (lx1/=size(f2,1) .or. lx2/=size(f2,2) .or. lx3/=size(f2,3) .or. lx1/=size(f3,1) .or. &
    lx2/=size(f3,2) .or. lx3/=size(f3,3)) error stop '!!! bad component sizes'


!ERROR CHECKING TO MAKE SURE DIFFRENCING IS DONE OVER A CONSISTENTLY-SIZED GRID
if (lx1 /= ubnd1-lbnd1+1 .or. lx2 /= ubnd2-lbnd2+1 .or. lx3 /= ubnd3-lbnd3+1) then
  error stop '!!!  Inconsistent array and mesh sizes in div3D function.'   !just bail on it and let the user figure it out
end if


!CHOOSE THE METRIC FACTORS VARIABLES BASED ON THE SIZE OF THE X3-VARIABLE, ALSO RECAST SO THE
!INDICES USED FOR F CAN ALSO BE USED IN THE METRIC FACTOR AND DX VARIABLE
!Can avoid wasting memory and copying of metric factor arrays by recoding with pointers
!Unfortunately a pointer is not guaranteed to be contiguous in memory so there may be a performance hit in doing this...
dx1=>x%dx1(lbnd1:ubnd1)
if (lx3<=x%lx3+4 .and. lx2<=x%lx2+4) then     !+4 in case we need to differentiate over ghost cells, e.g. in compression terms
  h1=>x%h1(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3(lbnd3:ubnd3)
  dx2=>x%dx2(lbnd2:ubnd2)
else if (lx3<=x%lx3all+4 .and. lx2<=x%lx2all+4) then     !presumes root or some process that has access to ALL full grid variables (normally only root).
  h1=>x%h1all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h2=>x%h2all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  h3=>x%h3all(lbnd1:ubnd1,lbnd2:ubnd2,lbnd3:ubnd3)
  dx3=>x%dx3all(lbnd3:ubnd3)
  dx2=>x%dx2all(lbnd2:ubnd2)
  print *,   '! Accessing root-only grid information in divergence function div3D'
else
  error stop '!!!  Array size is larger than full mesh or inconsistent sizes in data to be differentiated.'
end if

!REASSIGNING POINTERS AS SUBARRAYS TO BE DIFFERENCED PRODUCED CLEAN CODE (BELOW), BUT NOW REQUIRES THAT WE OFFSET
!ARRAY INDICES FOR THE POINTER VARIABLES, SINCE, BY DEFAULT, THEY START AT 1...  ACTUALLY THIS SHOULD NOT BE THE CASE
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

if (x%lx2>1) then   !only if the x2-direction is not null, this should be based on the grid data and not the data passed into this function
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


if (x%lx3>1) then    !only if non-singleton 3rd dimension (I'm not sure this will ever be the case due to dimension swapping)
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

end procedure div3D_curv_23

end submodule div
