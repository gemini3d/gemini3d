module advec

use phys_consts, only: lsp,ms

implicit none (external)

contains


subroutine advec_prep(isp,ns,rhovs1,vs1,vs2,vs3,rhoes,v1i,v2i,v3i)

!------------------------------------------------------------
!-------COMPUTE INTERFACE VELOCITIES AND LOAD UP GHOST CELLS
!------------------------------------------------------------
!-------Note that it is done on a per species basis

integer, intent(in) :: isp
real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,vs1,vs2,vs3,rhoes

real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4), intent(out) :: v1i   !why is size specified here???
real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4), intent(out) :: v2i
real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3), intent(out) :: v3i

real(wp), parameter :: vellim=2000.0
!    real(wp), parameter :: vellim=0.0
real(wp) :: coeff
integer :: ix2,ix3,lx1,lx2,lx3

lx1=size(vs1,1)-4
lx2=size(vs1,2)-4
lx3=size(vs1,3)-4


!COMPUTE INTERFACE VELCOTIES AND APPLY LIMITING, IF NEEDED
v1i(1,:,:)=vs1(1,1:lx2,1:lx3,isp)
v1i(2:lx1,:,:)=0.5*(vs1(1:lx1-1,1:lx2,1:lx3,isp)+vs1(2:lx1,1:lx2,1:lx3,isp))
!    v1i(lx1+1,:,:)=v1i(lx1,:,:)    !avoids issues with top boundary velocity spikes which may arise
v1i(lx1+1,:,:)=max(v1i(lx1,1:lx2,1:lx3),0.0)
v2i(:,1,:)=vs2(1:lx1,1,1:lx3,isp)
v2i(:,2:lx2,:)=0.5*(vs2(1:lx1,1:lx2-1,1:lx3,isp)+vs2(1:lx1,2:lx2,1:lx3,isp))
v2i(:,lx2+1,:)=vs2(1:lx1,lx2,1:lx3,isp)
v3i(:,:,1)=vs3(1:lx1,1:lx2,1,isp)
v3i(:,:,2:lx3)=0.5*(vs3(1:lx1,1:lx2,1:lx3-1,isp)+vs3(1:lx1,1:lx2,2:lx3,isp))
v3i(:,:,lx3+1)=vs3(1:lx1,1:lx2,lx3,isp)

!    if (isp<lsp-1) then
!      v1i(lx1+1,:,:)=max(vs1(lx1+1,:,:,isp),-1*vellim)    !limit inflow
!    else if (isp==6) then
!      v1i(lx1+1,:,:)=min(vs1(lx1+1,:,:,isp),10.0*vellim)      !limit outflow
!      v1i(lx1+1,:,:)=max(v1i(lx1+1,:,:),0.0)                 !limit inflow
!!        vadvx1(1,:,6)=mean(vadvx1(1,:,6));                          !for eq sims
!    else
!      v1i(lx1+1,:,:)=2.0*v1i(lx1,:,:)-v1i(lx1-1,:,:);      !cleans up large current situations
!    end if


!GHOST CELL VALUES FOR DENSITY
ns(:,0,:,isp)=ns(:,1,:,isp)
ns(:,-1,:,isp)=ns(:,1,:,isp)
ns(:,lx2+1,:,isp)=ns(:,lx2,:,isp)
ns(:,lx2+2,:,isp)=ns(:,lx2,:,isp)

ns(:,:,0,isp)=ns(:,:,1,isp)
ns(:,:,-1,isp)=ns(:,:,1,isp)
ns(:,:,lx3+1,isp)=ns(:,:,lx3,isp)
ns(:,:,lx3+2,isp)=ns(:,:,lx3,isp)

do ix3=1,lx3
  do ix2=1,lx2
    !logical bottom
    coeff=ns(2,ix2,ix3,isp)/ns(3,ix2,ix3,isp)
    ns(0,ix2,ix3,isp)=min(coeff*ns(1,ix2,ix3,isp),ns(1,ix2,ix3,isp))
    ns(-1,ix2,ix3,isp)=min(coeff*ns(0,ix2,ix3,isp),ns(0,ix2,ix3,isp))

    !logical top
    coeff=ns(lx1-1,ix2,ix3,isp)/ns(lx1-2,ix2,ix3,isp)
    ns(lx1+1,ix2,ix3,isp)=min(coeff*ns(lx1,ix2,ix3,isp),ns(lx1,ix2,ix3,isp))
    ns(lx1+2,ix2,ix3,isp)=min(coeff*ns(lx1+1,ix2,ix3,isp),ns(lx1+1,ix2,ix3,isp))
  end do
end do


!FOR X1 MOMENTUM DENSITY
rhovs1(:,0,:,isp)=rhovs1(:,1,:,isp);
rhovs1(:,-1,:,isp)=rhovs1(:,1,:,isp);
rhovs1(:,lx2+1,:,isp)=rhovs1(:,lx2,:,isp);
rhovs1(:,lx2+2,:,isp)=rhovs1(:,lx2,:,isp);

rhovs1(:,:,0,isp)=rhovs1(:,:,1,isp);
rhovs1(:,:,-1,isp)=rhovs1(:,:,1,isp);
rhovs1(:,:,lx3+1,isp)=rhovs1(:,:,lx3,isp);
rhovs1(:,:,lx3+2,isp)=rhovs1(:,:,lx3,isp);

rhovs1(0,:,:,isp)=2.0*v1i(1,:,:)-vs1(1,:,:,isp);  !initially these are velocities.  Also loose definition of 'conformable'
rhovs1(-1,:,:,isp)=rhovs1(0,:,:,isp)+rhovs1(0,:,:,isp)-vs1(1,:,:,isp);
rhovs1(lx1+1,:,:,isp)=2.0*v1i(lx1+1,:,:)-vs1(lx1,:,:,isp);
rhovs1(lx1+2,:,:,isp)=rhovs1(lx1+1,:,:,isp)+rhovs1(lx1+1,:,:,isp)-vs1(lx1,:,:,isp);

rhovs1(-1:0,:,:,isp)=rhovs1(-1:0,:,:,isp)*ns(-1:0,:,:,isp)*ms(isp)   !now convert to momentum density
rhovs1(lx1+1:lx1+2,:,:,isp)=rhovs1(lx1+1:lx1+2,:,:,isp)*ns(lx1+1:lx1+2,:,:,isp)*ms(isp)


!FOR INTERNAL ENERGY
rhoes(:,0,:,isp)=rhoes(:,1,:,isp);
rhoes(:,-1,:,isp)=rhoes(:,1,:,isp);
rhoes(:,lx2+1,:,isp)=rhoes(:,lx2,:,isp);
rhoes(:,lx2+2,:,isp)=rhoes(:,lx2,:,isp);

rhoes(:,:,0,isp)=rhoes(:,:,1,isp);
rhoes(:,:,-1,isp)=rhoes(:,:,1,isp);
rhoes(:,:,lx3+1,isp)=rhoes(:,:,lx3,isp);
rhoes(:,:,lx3+2,isp)=rhoes(:,:,lx3,isp);

rhoes(0,:,:,isp)=rhoes(1,:,:,isp);
rhoes(-1,:,:,isp)=rhoes(1,:,:,isp);
rhoes(lx1+1,:,:,isp)=rhoes(lx1,:,:,isp);
rhoes(lx1+2,:,:,isp)=rhoes(lx1,:,:,isp);

end subroutine advec_prep


function advec3D_MC(f,v1i,v2i,v3i,dt,dx1,dx1i,dx2,dx2i,dx3,dx3i)
real(wp), dimension(-1:,-1:,-1:), intent(in) :: f    !f includes ghost cells
real(wp), dimension(:,:,:), intent(in) :: v1i
real(wp), dimension(0:), intent(in) :: dx1       !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i
real(wp), dimension(:,:,:), intent(in) :: v2i
real(wp), dimension(0:), intent(in) :: dx2
real(wp), dimension(:), intent(in) :: dx2i
real(wp), dimension(:,:,:), intent(in) :: v3i
real(wp), dimension(0:), intent(in) :: dx3
real(wp), dimension(:), intent(in) :: dx3i
real(wp), intent(in) :: dt

integer :: ix1,ix2,ix3,lx1,lx2,lx3
real(wp), dimension(-1:size(f,1)-2) :: fx1slice
real(wp), dimension(1:size(f,1)-3) :: v1slice
real(wp), dimension(-1:size(f,2)-2) :: fx2slice
real(wp), dimension(1:size(f,2)-3) :: v2slice
real(wp), dimension(-1:size(f,3)-2) :: fx3slice
real(wp), dimension(1:size(f,3)-3) :: v3slice
real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: advec3D_MC

lx1=size(f,1)-4
lx2=size(f,2)-4
lx3=size(f,3)-4

!x1-sweep
do ix3=1,lx3
  do ix2=1,lx2
    fx1slice=f(:,ix2,ix3);
    v1slice=v1i(:,ix2,ix3);
    fx1slice=advec1D_MC(fx1slice,v1slice,dt,dx1,dx1i)
    advec3D_MC(:,ix2,ix3)=fx1slice;
  end do
end do

!copy x2,x3 boundary conditions to partially updated variable for next sweeps
advec3D_MC(:,-1:0,:)=f(:,-1:0,:);
advec3D_MC(:,lx2+1:lx2+2,:)=f(:,lx2+1:lx2+2,:);
advec3D_MC(:,:,-1:0)=f(:,:,-1:0);
advec3D_MC(:,:,lx3+1:lx3+2)=f(:,:,lx3+1:lx3+2);

!x2-sweep
do ix3=1,lx3
  do ix1=1,lx1
    fx2slice=advec3D_MC(ix1,:,ix3)
    v2slice=v2i(ix1,:,ix3)
    fx2slice=advec1D_MC(fx2slice,v2slice,dt,dx2,dx2i)
    advec3D_MC(ix1,:,ix3)=fx2slice
  end do
end do

!x3-sweep
do ix2=1,lx2
  do ix1=1,lx1
    fx3slice=advec3D_MC(ix1,ix2,:)
    v3slice=v3i(ix1,ix2,:)
    fx3slice=advec1D_MC(fx3slice,v3slice,dt,dx3,dx3i)
    advec3D_MC(ix1,ix2,:)=fx3slice
  end do
end do

end function advec3D_MC



function advec2D_MC(f,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
real(wp), dimension(-1:,-1:), intent(in) :: f    !f includes ghost cells
real(wp), dimension(:,:), intent(in) :: v1i
real(wp), dimension(0:), intent(in) :: dx1       !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i
real(wp), dimension(:,:), intent(in) :: v2i
real(wp), dimension(0:), intent(in) :: dx2
real(wp), dimension(:), intent(in) :: dx2i
real(wp), intent(in) :: dt

integer :: ix1,ix2,lx1,lx2
real(wp), dimension(-1:size(f,1)-2) :: fx1slice
real(wp), dimension(1:size(f,1)-3) :: v1slice
real(wp), dimension(-1:size(f,2)-2) :: fx2slice
real(wp), dimension(1:size(f,2)-3) :: v2slice
real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2) :: advec2D_MC

lx1=size(f,1)-4
lx2=size(f,2)-4

!x1-sweep
do ix2=1,lx2
  fx1slice=f(:,ix2);
  v1slice=v1i(:,ix2);
  fx1slice=advec1D_MC(fx1slice,v1slice,dt,dx1,dx1i)
  advec2D_MC(:,ix2)=fx1slice;
end do

!copy x2 boundary conditions to partially updated variable for next sweep
advec2D_MC(:,-1:0)=f(:,-1:0);
advec2D_MC(:,lx2+1:lx2+2)=f(:,lx2+1:lx2+2);

!x2-sweep (it appears that fortran doesn't differentiate b/t row and column arrays, so no reshapes...)
do ix1=1,lx1
  fx2slice=advec2D_MC(ix1,:)
  v2slice=v2i(ix1,:)
  fx2slice=advec1D_MC(fx2slice,v2slice,dt,dx2,dx2i)
  advec2D_MC(ix1,:)=fx2slice
end do
end function advec2D_MC



function advec1D_MC(f,v1i,dt,dx1,dx1i)
real(wp), dimension(-1:), intent(in) :: f    !f includes ghost cells
real(wp), dimension(:), intent(in) :: v1i
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i
real(wp), intent(in) :: dt

integer :: ix1,lx1
real(wp), dimension(size(v1i)) :: phi
real(wp), dimension(0:size(v1i)) :: slope         !slopes only need through first layer of ghosts
real(wp) :: lslope,rslope,cslope
real(wp), dimension(-1:size(f)-2) :: advec1D_MC


lx1=size(f)-4

!Slopes
lslope=(f(0)-f(-1))/dx1(0)
do ix1=0,lx1+1
  rslope=(f(ix1+1)-f(ix1))/dx1(ix1+1)
  cslope=(f(ix1+1)-f(ix1-1))/(dx1(ix1)+dx1(ix1+1))
  slope(ix1)=minmod(cslope,minmod(2*lslope,2*rslope))

  lslope=rslope
end do


!Slope-limited flux at ***left*** wall of cell ix1.
!The treatment of slope here (ie the dx1(ix1)) assumes that the grid point isp centered within the cell
do ix1=1,lx1+1
  if (v1i(ix1) < 0) then
    phi(ix1)=f(ix1)*v1i(ix1) - 0.5*v1i(ix1)*(dx1(ix1)+v1i(ix1)*dt)*slope(ix1)
  else
    phi(ix1)=f(ix1-1)*v1i(ix1) + 0.5*v1i(ix1)*(dx1(ix1)-v1i(ix1)*dt)*slope(ix1-1)
  end if
end do


!flux differencing form
advec1D_MC(1:lx1)=f(1:lx1)-dt*(phi(2:lx1+1)-phi(1:lx1))/dx1i


!default to copying ghost cells
advec1D_MC(-1:0)=f(-1:0)
advec1D_MC(lx1+1:lx1+2)=f(lx1+1:lx1+2)
end function advec1D_MC



elemental real(wp) function minmod(a,b)
real(wp), intent(in) :: a,b

if (a*b <= 0.) then
  minmod=0.
else if (abs(a) < abs(b)) then
  minmod=a
else
  minmod=b
end if
end function minmod



pure function advec1D_DC(f,v1i,dt,dx1,dx1i)
real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
real(wp), dimension(:), intent(in) :: dx1i
real(wp), dimension(-1:), intent(in) :: f
real(wp), dimension(:), intent(in) :: v1i
real(wp), intent(in) :: dt

integer :: ix1,lx1
real(wp), dimension(size(v1i)) :: phi
real(wp), dimension(-1:size(f)-2) :: advec1D_DC


lx1=size(f)-4 !size of main grid excluding ghost cells

do ix1=1,lx1+1
  if (v1i(ix1) < 0) then
    phi(ix1)=f(ix1)*v1i(ix1)
  else
    phi(ix1)=f(ix1-1)*v1i(ix1)
  end if
end do

advec1D_DC(1:lx1)=f(1:lx1)-dt*(phi(2:lx1+1)-phi(1:lx1))/dx1i
advec1D_DC(-1:0)=f(-1:0)
advec1D_DC(lx1+1:lx1+2)=f(lx1+1:lx1+2)
end function advec1D_DC

end module advec
