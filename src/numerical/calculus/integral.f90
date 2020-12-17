submodule (calculus) integral

implicit none (type, external)

contains

module procedure integral3D1_curv
!  integral3D1_curv(f,x,lbnd,ubnd)
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
integer :: ix1,lx1

lx1=size(f,1)

integral3D1_curv(1,:,:)=0._wp
do ix1=2,lx1
  integral3D1_curv(ix1,:,:)=integral3D1_curv(ix1-1,:,:)+0.5_wp*(f(ix1,:,:)+f(ix1-1,:,:))*x%dx1(lbnd+ix1-1)
end do
end procedure integral3D1_curv


module procedure integral3D1_curv_alt

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
integer :: ix1,lx1

lx1=size(f,1)

integral3D1_curv_alt(lx1,:,:)=0._wp
do ix1=lx1-1,1,-1
  !! start from the logical top and sum  downward (keep dx positive since the intent is to integrate from fixed point to top)
  integral3D1_curv_alt(ix1,:,:)=integral3D1_curv_alt(ix1+1,:,:)+0.5_wp*(f(ix1,:,:)+f(ix1+1,:,:))*x%dx1(lbnd+ix1-1+1)
  !! +1 since we are starting from top
end do
end procedure integral3D1_curv_alt


module procedure integral2D1_curv

!PURPOSEFULLY KEPT WARNING TO MAINTIN CONGRUITY WITH OTHER CALCULUS FUNCTIONS
!/home/zettergm/zettergmdata/GEMINI/numerical/calculus/calculus.f90:568:0:
!warning: unused parameter ‘ubnd’ [-Wunused-parameter]
!   function integral2D1_curv(f,x,lbnd,ubnd)

!------------------------------------------------------------
!-------COMPUTE AN INTEGRAL OF A 2D ARRAY ALONG THE 1-DIM.  IT IS EXPECTED THAT
!-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
!-------THEY ARE PASSED INTO THIS ROUTINE
!------------------------------------------------------------
integer :: ix1,lx1

lx1=size(f,1)

integral2D1_curv(1,:)=0._wp
do ix1=2,lx1
  integral2D1_curv(ix1,:)=integral2D1_curv(ix1-1,:)+0.5_wp*(f(ix1,:)+f(ix1-1,:))*x%dx1(lbnd+ix1-1)
end do
end procedure integral2D1_curv


module procedure integral2D1_curv_alt

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
integer :: ix2,lx2

lx2=size(f,1)

integral2D1_curv_alt(1,:)=0._wp
do ix2=2,lx2
  integral2D1_curv_alt(ix2,:)=integral2D1_curv_alt(ix2-1,:)+0.5_wp*(f(ix2,:)+f(ix2-1,:))*x%dx2(lbnd+ix2-1)
end do
end procedure integral2D1_curv_alt


module procedure integral2D2_curv

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

integer :: ix2,lx2

lx2=size(f,2)

integral2D2_curv(:,1)=0._wp
do ix2=2,lx2
  integral2D2_curv(:,ix2)=integral2D2_curv(:,ix2-1)+0.5_wp*(f(:,ix2)+f(:,ix2-1))*x%dx2all(lbnd+ix2-1)
end do
end procedure integral2D2_curv


module procedure integral2D2_curv_alt

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
integer :: ix3,lx3

lx3=size(f,2)

integral2D2_curv_alt(:,1)=0._wp
do ix3=2,lx3
  integral2D2_curv_alt(:,ix3)=integral2D2_curv_alt(:,ix3-1)+0.5_wp*(f(:,ix3)+f(:,ix3-1))*x%dx3all(lbnd+ix3-1)
end do
end procedure integral2D2_curv_alt

end submodule integral
