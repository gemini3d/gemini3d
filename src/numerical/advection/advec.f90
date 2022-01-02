module advec

!> this module contains advection-related procedures that are independent on mpi library calls

use phys_consts, only: lsp,ms,wp
use grid, only : gridflag
use meshobj, only: curvmesh
  !! do not import grid sizes in case we want do subgrid advection...

implicit none (type, external)
private
public :: interface_vels_allspec,sweep3_allspec,sweep1_allspec,sweep2_allspec

contains
  !> do averaging to compute cell interface velocities; requires pre-haloing in order for cells
  !    to have updated boundary info from ghost cells.  Single species only, "isp"; does not require
  !    any knowledge about the mpi image configuration(s).
  subroutine interface_vels(isp,vs1,vs2,vs3,v1i,v2i,v3i)
    integer, intent(in) :: isp
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4), intent(inout) :: v1i
    !! intent(out)
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4), intent(inout) :: v2i
    !! intent(out)
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3), intent(inout) :: v3i
    integer :: lx1,lx2,lx3
  
    lx1=size(vs1,1)-4
    lx2=size(vs1,2)-4
    lx3=size(vs1,3)-4
    
    !COMPUTE INTERFACE VELCOTIES AND APPLY LIMITING, IF NEEDED
    v1i(2:lx1,:,:)=0.5*(vs1(1:lx1-1,1:lx2,1:lx3,isp)+vs1(2:lx1,1:lx2,1:lx3,isp))   !first the interior points
    if (gridflag==0) then          !closed dipole grid
      v1i(1,:,:)=vs1(1,1:lx2,1:lx3,isp)
      v1i(lx1+1,:,:)=vs1(lx1,1:lx2,1:lx3,isp)
    else if (gridflag==1) then     !inverted grid (assumes northern hemisphere???)
      v1i(lx1+1,:,:)=vs1(lx1,1:lx2,1:lx3,isp)
      !! lowest alt on grid.
      v1i(1,:,:) = min(v1i(2,1:lx2,1:lx3), 0._wp)
      !! highest alt; interesting that this is not vs1...
    else                           !some type of non-inverted grid
      v1i(1,:,:) = vs1(1,1:lx2,1:lx3,isp)
    !!    v1i(lx1+1,:,:)=v1i(lx1,:,:)    !avoids issues with top boundary velocity spikes which may arise
      v1i(lx1+1,:,:) = max(v1i(lx1,1:lx2,1:lx3),0._wp)
      !! NOTE: interesting that this is not vs1...
    end if
  
    ! THIS TYPE OF LIMITING MAY BE NEEDED FOR VERY HIGH-ALTITUDE SIMULATIONS...
    !    if (isp<lsp-1) then
    !      v1i(lx1+1,:,:)=max(vs1(lx1+1,:,:,isp),-1*vellim)    !limit inflow
    !    else if (isp==6) then
    !      v1i(lx1+1,:,:)=min(vs1(lx1+1,:,:,isp),10.0*vellim)      !limit outflow
    !      v1i(lx1+1,:,:)=max(v1i(lx1+1,:,:),0.0)                 !limit inflow
    !!        vadvx1(1,:,6)=mean(vadvx1(1,:,6))                          !for eq sims
    !    else
    !      v1i(lx1+1,:,:)=2.0*v1i(lx1,:,:)-v1i(lx1-1,:,:)      !cleans up large current situations
    !    end if
  
    !> AFTER HALOING CAN COMPUTE THE X3 INTERFACE VELOCITIES NORMALLY
    v2i(:,1:lx2+1,:)=0.5_wp*(vs2(1:lx1,0:lx2,1:lx3,isp)+vs2(1:lx1,1:lx2+1,1:lx3,isp))
    v3i(:,:,1:lx3+1)=0.5_wp*(vs3(1:lx1,1:lx2,0:lx3,isp)+vs3(1:lx1,1:lx2,1:lx3+1,isp))
  end subroutine interface_vels
  
  
  !> compute cell interface velocities for all species being simulated
  subroutine interface_vels_allspec(vs1,vs2,vs3,vs1i,vs2i,vs3i,lsp)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4,1:size(vs1,4)), intent(inout) :: vs1i
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4,1:size(vs2,4)), intent(inout) :: vs2i
    real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3,1:size(vs3,4)), intent(inout) :: vs3i
    integer, intent(in) :: lsp
    integer :: isp
  
    if (lsp>size(vs1,4)) error stop 'number of interface vels must be less than or equal to total species number'
    do isp=1,lsp
      call interface_vels(isp,vs1,vs2,vs3,vs1i(:,:,:,isp),vs2i(:,:,:,isp),vs3i(:,:,:,isp))
    end do
  end subroutine interface_vels_allspec


  !> sweep all species along the 1 axis
  subroutine sweep1_allspec(fs,vs1i,dt,x,lsp)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: fs    !fs includes ghost cells and all species
    real(wp), dimension(:,:,:,:), intent(in) :: vs1i           ! includes all species velocities
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: lsp                                 ! sweep the first "lsp" species only
    integer :: isp
  
    if (lsp>size(fs,4)) error stop 'number of swept species must be less than or equal to total species number'
    do isp=1,lsp
      call sweep1(fs(:,:,:,isp),vs1i(:,:,:,isp),dt,x)
    end do
  end subroutine sweep1_allspec
  
  
  !> 2-dimensionally split transport for all species
  subroutine sweep2_allspec(fs,vs2i,dt,x,frank,lsp)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: fs    !fs includes ghost cells and all species
    real(wp), dimension(:,:,:,:), intent(in) :: vs2i           ! includes all species velocities
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: frank
    integer, intent(in) :: lsp
    integer :: isp
  
    if (lsp>size(fs,4)) error stop 'number of swept species must be less than or equal to total species number'
    do isp=1,lsp
      call sweep2(fs(:,:,:,isp),vs2i(:,:,:,isp),dt,x,frank)
    end do
  end subroutine sweep2_allspec
  
  
  !> 3-dimensionally split transport for all species
  subroutine sweep3_allspec(fs,vs3i,dt,x,frank,lsp)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: fs    !fs includes ghost cells and all species
    real(wp), dimension(:,:,:,:), intent(in) :: vs3i           ! includes all species velocities
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: frank
    integer, intent(in) :: lsp
    integer :: isp
  
    if (lsp>size(fs,4)) error stop 'number of swept species must be less than or equal to total species number'
    do isp=1,lsp
      call sweep3(fs(:,:,:,isp),vs3i(:,:,:,isp),dt,x,frank)
    end do
  end subroutine sweep3_allspec


  !> dimensionally split advection along the 1-axis for a single species (i.e. only 3D arrays)
  subroutine sweep1(f,v1i,dt,x)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: f    !f includes ghost cells
    real(wp), dimension(:,:,:), intent(in) :: v1i
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    real(wp), dimension(-1:size(f,1)-2) :: fx1slice
    real(wp), dimension(1:size(f,1)-3) :: v1slice
    real(wp), dimension(-1:size(f,1)-2) :: h11x1slice    !includes ghost cells
    real(wp), dimension(1:size(f,1)-3) :: h12ix1slice    !just includes interface info
    real(wp), dimension(1:size(f,1)-3) :: h1ix1slice
    integer :: ix2,ix3,lx2,lx3
  
    lx2=size(f,2)-4
    lx3=size(f,3)-4
  
    do ix3=1,lx3
      do ix2=1,lx2
        fx1slice=f(:,ix2,ix3)
        v1slice=v1i(:,ix2,ix3)
        h11x1slice=x%h1(:,ix2,ix3)*x%h2(:,ix2,ix3)*x%h3(:,ix2,ix3)
        h12ix1slice=x%h2x1i(:,ix2,ix3)*x%h3x1i(:,ix2,ix3)
        h1ix1slice=x%h1x1i(:,ix2,ix3)
        fx1slice=advec1D_MC_curv(fx1slice,v1slice,dt,x%dx1,x%dx1i,h11x1slice,h12ix1slice,h1ix1slice)
        f(:,ix2,ix3)=fx1slice;
      end do
    end do
  end subroutine sweep1
  
  
  !> Dimensionally split advection in the second direction/axis
  subroutine sweep2(f,v2i,dt,x,frank)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: f    !f includes ghost cells
    real(wp), dimension(:,:,:), intent(in) :: v2i
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: frank
    real(wp), dimension(-1:size(f,2)-2) :: fx2slice
    real(wp), dimension(1:size(f,2)-3) :: v2slice
    real(wp), dimension(-1:size(f,2)-2) :: h21x2slice    !includes ghost cells
    real(wp), dimension(1:size(f,2)-3) :: h22ix2slice    !just includes interface info
    real(wp), dimension(1:size(f,2)-3) :: h2ix2slice
    integer :: ix1,ix3,lx1,lx3
  
    lx1=size(f,1)-4
    lx3=size(f,3)-4
  
    do ix3=1,lx3
      do ix1=1,lx1
        fx2slice=f(ix1,:,ix3)
        v2slice=v2i(ix1,:,ix3)
        if (frank==0) then
          h21x2slice=x%h1(ix1,:,ix3)*x%h2(ix1,:,ix3)*x%h3(ix1,:,ix3)
          h22ix2slice=x%h1x2i(ix1,:,ix3)*x%h3x2i(ix1,:,ix3)
        else
          h21x2slice=x%h1(ix1,:,ix3)**2*x%h2(ix1,:,ix3)*x%h3(ix1,:,ix3)
          h22ix2slice=x%h1x2i(ix1,:,ix3)**2*x%h3x2i(ix1,:,ix3)
        end if
        h2ix2slice=x%h2x2i(ix1,:,ix3)
        fx2slice=advec1D_MC_curv(fx2slice,v2slice,dt,x%dx2,x%dx2i,h21x2slice,h22ix2slice,h2ix2slice)
        f(ix1,:,ix3)=fx2slice
      end do
    end do
  end subroutine sweep2
  
  
  !> Do a dimensionally split advection sweep along the third dimension/axis
  subroutine sweep3(f,v3i,dt,x,frank)
    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: f    !f includes ghost cells
    real(wp), dimension(:,:,:), intent(in) :: v3i
    real(wp), intent(in) :: dt
    class(curvmesh), intent(in) :: x
    integer, intent(in) :: frank
    real(wp), dimension(-1:size(f,3)-2) :: fx3slice
    real(wp), dimension(1:size(f,3)-3) :: v3slice
    real(wp), dimension(-1:size(f,3)-2) :: h31x3slice    !includes ghost cells
    real(wp), dimension(1:size(f,3)-3) :: h32ix3slice    !just includes interface info
    real(wp), dimension(1:size(f,3)-3) :: h3ix3slice
    integer :: ix1,ix2,lx1,lx2
  
    lx1=size(f,1)-4
    lx2=size(f,2)-4
  
    do ix2=1,lx2
      do ix1=1,lx1
        fx3slice=f(ix1,ix2,:)
        v3slice=v3i(ix1,ix2,:)
        if (frank==0) then     !advecting a scalar quantity
          h31x3slice=x%h1(ix1,ix2,:)*x%h2(ix1,ix2,:)*x%h3(ix1,ix2,:)
          h32ix3slice=x%h1x3i(ix1,ix2,:)*x%h2x3i(ix1,ix2,:)
        else    !advecting 1-component of a rank 1 tensor
          h31x3slice=x%h1(ix1,ix2,:)**2*x%h2(ix1,ix2,:)*x%h3(ix1,ix2,:)
          h32ix3slice=x%h1x3i(ix1,ix2,:)**2*x%h2x3i(ix1,ix2,:)
        end if
        h3ix3slice=x%h3x3i(ix1,ix2,:)
        fx3slice=advec1D_MC_curv(fx3slice,v3slice,dt,x%dx3,x%dx3i,h31x3slice,h32ix3slice,h3ix3slice)
        f(ix1,ix2,:)=fx3slice
      end do
    end do
  end subroutine sweep3


  function advec1D_MC_curv(f,v1i,dt,dx1,dx1i,ha1,ha2i,h1i)
  !----------------------------------------------------------------
  !---- Generic advection routine, can account for metric factors
  !-----  via arguments including for vector quantities.
  !----------------------------------------------------------------
    real(wp), dimension(-1:), intent(in) :: f    !f includes ghost cells
    real(wp), dimension(:), intent(in) :: v1i
    real(wp), intent(in) :: dt
    real(wp), dimension(0:), intent(in) :: dx1   !ith backward difference
    real(wp), dimension(:), intent(in) :: dx1i   !interface-based differences - does not include any ghost cell
    real(wp), dimension(-1:), intent(in) :: ha1    !cell-centered metric factor product 1; includes ghost cells
    real(wp), dimension(:), intent(in) :: ha2i   !cell interface metric factor product 2
    real(wp), dimension(:), intent(in) :: h1i    !cell interface metric factor for dimension being advected
    integer :: ix1,lx1                          !overwrite grid module lx1 in this function's scope, since it gets called for x1,x2,x3
    real(wp), dimension(size(v1i)) :: phi
    real(wp), dimension(0:size(v1i)) :: slope         !slopes only need through first layer of ghosts
    real(wp) :: lslope,rslope,cslope
    real(wp), dimension(-1:size(f)-2) :: advec1D_MC_curv
    
    lx1=size(f)-4     ! we don't know what dimension this is so we actually do need to compute the size
  
    if (lx1>1) then     ! don't advect a single computational point
      !Slopes
      lslope=(f(0)-f(-1))/dx1(0)
      do ix1=0,lx1+1
        rslope=(f(ix1+1)-f(ix1))/dx1(ix1+1)
        cslope=(f(ix1+1)-f(ix1-1))/(dx1(ix1)+dx1(ix1+1))
        slope(ix1)=minmod(cslope,minmod(2*lslope,2*rslope))
      
        lslope=rslope
      end do
      
      !Slope-limited flux at ***left*** wall of cell ix1.
      !The treatment of slope here (ie the dx1(ix1)) assumes that the grid point is centered within the cell
      do ix1=1,lx1+1
        if (v1i(ix1) < 0) then
          phi(ix1)=f(ix1)*v1i(ix1) - 0.5_wp*v1i(ix1)*(dx1(ix1)+v1i(ix1)/h1i(ix1)*dt)*slope(ix1)
        else
          phi(ix1)=f(ix1-1)*v1i(ix1) + 0.5_wp*v1i(ix1)*(dx1(ix1)-v1i(ix1)/h1i(ix1)*dt)*slope(ix1-1)
        end if
      end do
      
      !flux differencing form
      advec1D_MC_curv(1:lx1)=f(1:lx1)-dt*(ha2i(2:lx1+1)*phi(2:lx1+1)-ha2i(1:lx1)*phi(1:lx1))/dx1i/ha1(1:lx1)
      
      !default to copying ghost cells
      advec1D_MC_curv(-1:0)=f(-1:0)
      advec1D_MC_curv(lx1+1:lx1+2)=f(lx1+1:lx1+2)
    else
      advec1D_MC_curv(:)=f(:)
    end if
  end function advec1D_MC_curv
  
  
  elemental real(wp) function minmod(a,b)
    real(wp), intent(in) :: a,b
    if (a*b <= 0._wp) then
      minmod = 0._wp
    else if (abs(a) < abs(b)) then
      minmod=a
    else
      minmod=b
    end if
  end function minmod
end module advec
