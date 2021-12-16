module advec_mpi

use phys_consts, only: lsp,ms, wp
use grid, only : gridflag
use meshobj, only: curvmesh
  !! do not import grid sizes in case we want do subgrid advection...
use mpimod, only: mpi_cfg, halo, tag=>gemini_mpi

implicit none (type, external)
private
public :: advec3d_mc_mpi, advec_prep_mpi, halo_interface_vels_allspec,interface_vels_allspec,set_global_boundaries_allspec, &
          sweep3_allspec,sweep1_allspec,sweep2_allspec

!> OVERLOAD ADVECTION TO DEAL WITH THE CURVILINEAR GRID/MESH STRUCTURE.
!> NOTE THAT THE LOWER-LEVEL CALLS ARE DISTINCT, NOT-OVERLOADED PROCEDURES.
interface advec3D_MC_mpi
  module procedure advec3D_MC_mpi_curv_23
end interface advec3D_MC_mpi
interface advec_prep_mpi
  module procedure advec_prep_mpi_23
end interface advec_prep_mpi

contains


subroutine advec_prep_mpi_23(isp,isperiodic,ns,rhovs1,vs1,vs2,vs3,rhoes,v1i,v2i,v3i)
!! COMPUTE INTERFACE VELOCITIES AND LOAD UP GHOST CELLS
!! FOR FLUID STATE VARIABLES
!!
!! Note that it is done on a per species basis
!! 5/23/2015 - may need to be changed for 2D/1D sims which
!! have only one element in the x2 direction...
  integer, intent(in) :: isp
  logical, intent(in) :: isperiodic
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,vs1,vs2,vs3,rhoes
  real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4), intent(inout) :: v1i
  !! intent(out)
  real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4), intent(inout) :: v2i
  !! intent(out)
  real(wp), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3), intent(inout) :: v3i
  !! intent(out)
  !    real(wp), parameter :: vellim=2000.0
  !    real(wp), parameter :: vellim=0.0
  real(wp) :: coeff
  integer :: ix2,ix3,lx1,lx2,lx3
  integer :: idleft,idright,idup,iddown
  real(wp), dimension(-1:size(vs3,1)-2,-1:size(vs3,2)-2,-1:size(vs3,3)-2) :: param,param2,param3,param4
  real(wp) :: tstart,tfini
  
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
  
  !v2i(:,1,:)=vs2(1:lx1,1,1:lx3,isp)
  !v2i(:,2:lx2,:)=0.5*(vs2(1:lx1,1:lx2-1,1:lx3,isp)+vs2(1:lx1,2:lx2,1:lx3,isp))
  !v2i(:,lx2+1,:)=vs2(1:lx1,lx2,1:lx3,isp)
  
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
  
  !> GHOST CELL VALUES FOR DENSITY (these need to be done last to avoid being overwritten by send/recv's!!!)
  do ix3=1,lx3
    do ix2=1,lx2
      !> logical bottom
      coeff=ns(2,ix2,ix3,isp)/ns(3,ix2,ix3,isp)
      ns(0,ix2,ix3,isp)=min(coeff*ns(1,ix2,ix3,isp),ns(1,ix2,ix3,isp))
      ns(-1,ix2,ix3,isp)=min(coeff*ns(0,ix2,ix3,isp),ns(0,ix2,ix3,isp))
  
      !> logical top
      coeff=ns(lx1-1,ix2,ix3,isp)/ns(lx1-2,ix2,ix3,isp)
      ns(lx1+1,ix2,ix3,isp)=min(coeff*ns(lx1,ix2,ix3,isp),ns(lx1,ix2,ix3,isp))
      ns(lx1+2,ix2,ix3,isp)=min(coeff*ns(lx1+1,ix2,ix3,isp),ns(lx1+1,ix2,ix3,isp))
    end do
  end do
  
  !FOR X1 MOMENTUM DENSITY
  rhovs1(0,1:lx2,1:lx3,isp)=2*v1i(1,:,:)-vs1(1,1:lx2,1:lx3,isp)
  !! initially these are velocities.  Also loose definition of 'conformable'.  Also corners never get set, but I suppose they aren't really used anyway.
  rhovs1(-1,:,:,isp)=rhovs1(0,:,:,isp)+rhovs1(0,:,:,isp)-vs1(1,:,:,isp)
  rhovs1(lx1+1,1:lx2,1:lx3,isp)=2*v1i(lx1+1,:,:)-vs1(lx1,1:lx2,1:lx3,isp)
  rhovs1(lx1+2,:,:,isp)=rhovs1(lx1+1,:,:,isp)+rhovs1(lx1+1,:,:,isp)-vs1(lx1,:,:,isp)
  
  rhovs1(-1:0,:,:,isp)=rhovs1(-1:0,:,:,isp)*ns(-1:0,:,:,isp)*ms(isp)
  !! now convert to momentum density
  rhovs1(lx1+1:lx1+2,:,:,isp)=rhovs1(lx1+1:lx1+2,:,:,isp)*ns(lx1+1:lx1+2,:,:,isp)*ms(isp)
  
  !> FOR INTERNAL ENERGY
  rhoes(0,:,:,isp)=rhoes(1,:,:,isp)
  rhoes(-1,:,:,isp)=rhoes(1,:,:,isp)
  rhoes(lx1+1,:,:,isp)=rhoes(lx1,:,:,isp)
  rhoes(lx1+2,:,:,isp)=rhoes(lx1,:,:,isp)
  
  !MZ - collect the x2 boundary conditions here - these are no longer global
  iddown=mpi_cfg%myid2-1
  idup=mpi_cfg%myid2+1
  
  !> NEED TO ALSO PASS THE X2 VELOCITIES SO WE CAN COMPUTE INTERFACE VALUES
  param=vs2(:,:,:,isp)
  call halo(param,1,tag%vs2BC,isperiodic)
  !! we only need one ghost cell to compute interface velocities
  vs2(:,:,:,isp)=param
  
  !> SET THE GLOBAL X2 BOUNDARY CELLS AND ASSUME HALOING WON'T OVERWRITE.
  !> THIS DIMENSION IS ASSUMED TO NEVER BE PEREIODIC
  if (iddown==-1) then
    vs2(:,0,:,isp)=vs2(:,1,:,isp)
  
    ns(:,0,:,isp)=ns(:,1,:,isp)
    ns(:,-1,:,isp)=ns(:,1,:,isp)
    rhovs1(:,0,:,isp)=rhovs1(:,1,:,isp)
    rhovs1(:,-1,:,isp)=rhovs1(:,1,:,isp)
    rhoes(:,0,:,isp)=rhoes(:,1,:,isp)
    rhoes(:,-1,:,isp)=rhoes(:,1,:,isp)
  end if
  if (idup==mpi_cfg%lid2) then
    vs2(:,lx2+1,:,isp)=vs2(:,lx2,:,isp)
  
    ns(:,lx2+1,:,isp)=ns(:,lx2,:,isp)
    ns(:,lx2+2,:,isp)=ns(:,lx2,:,isp)
    rhovs1(:,lx2+1,:,isp)=rhovs1(:,lx2,:,isp)
    rhovs1(:,lx2+2,:,isp)=rhovs1(:,lx2,:,isp)
    rhoes(:,lx2+1,:,isp)=rhoes(:,lx2,:,isp)
    rhoes(:,lx2+2,:,isp)=rhoes(:,lx2,:,isp)
  end if
  
  !> NOW DEAL WITH ADVECTION ALONG X3; FIRST IDENTIFY MY NEIGHBORS
  idleft=mpi_cfg%myid3-1; idright=mpi_cfg%myid3+1
  
  !> PASS X3 VELOCITY BOUNDARY CONDITIONS WITH GENERIC HALOING ROUTINES
  param=vs3(:,:,:,isp)
  call halo(param,1,tag%vs3BC,isperiodic)
  !! we only need one ghost cell to compute interface velocities
  vs3(:,:,:,isp)=param
  
  !> these will now be haloed internal to the advection routines, viz. all advected quantities are haloed withine advection
  !param2=ns(:,:,:,isp)
  !call halo(param2,2,tag%nsBC)
  !ns(:,:,:,isp)=param2
  !
  !param3=rhovs1(:,:,:,isp)
  !call halo(param3,2,tag%rhovs1BC)
  !rhovs1(:,:,:,isp)=param3
  !
  !param4=rhoes(:,:,:,isp)
  !call halo(param4,2,tag%rhoesBC)
  !rhoes(:,:,:,isp)=param4
  
  !> SET THE GLOBAL x3 BOUNDARY CELLS AND ASSUME THAT HALOING WON'T OVERWRITE...
  if (.not. isperiodic) then
    if (idleft==-1) then
      !! left side is at global boundary, assume haloing won't overwrite
      vs3(:,:,0,isp)=vs3(:,:,1,isp)
      !! copy first cell to first ghost (vs3 not advected so only need only ghost)
  
      ns(:,:,0,isp)=ns(:,:,1,isp)
      ns(:,:,-1,isp)=ns(:,:,1,isp)
      rhovs1(:,:,0,isp)=rhovs1(:,:,1,isp)
      rhovs1(:,:,-1,isp)=rhovs1(:,:,1,isp)
      rhoes(:,:,0,isp)=rhoes(:,:,1,isp)
      rhoes(:,:,-1,isp)=rhoes(:,:,1,isp)
    end if
    if (idright==mpi_cfg%lid3) then    !my right boundary is the global boundary, assume haloing won't overwrite
      vs3(:,:,lx3+1,isp)=vs3(:,:,lx3,isp)    !copy last cell to first ghost (all that's needed since vs3 not advected)
  
      ns(:,:,lx3+1,isp)=ns(:,:,lx3,isp)
      ns(:,:,lx3+2,isp)=ns(:,:,lx3,isp)
      rhovs1(:,:,lx3+1,isp)=rhovs1(:,:,lx3,isp)
      rhovs1(:,:,lx3+2,isp)=rhovs1(:,:,lx3,isp)
      rhoes(:,:,lx3+1,isp)=rhoes(:,:,lx3,isp)
      rhoes(:,:,lx3+2,isp)=rhoes(:,:,lx3,isp)
    end if
  end if
  
  !> AFTER HALOING CAN COMPUTE THE X3 INTERFACE VELOCITIES NORMALLY
  v2i(:,1:lx2+1,:)=0.5_wp*(vs2(1:lx1,0:lx2,1:lx3,isp)+vs2(1:lx1,1:lx2+1,1:lx3,isp))
  v3i(:,:,1:lx3+1)=0.5_wp*(vs3(1:lx1,1:lx2,0:lx3,isp)+vs3(1:lx1,1:lx2,1:lx3+1,isp))
end subroutine advec_prep_mpi_23


!> Perform haloing needed to ghost-fill so cell interface vels (single species) can be computed across the grid
subroutine halo_interface_vels(isp,isperiodic,vs2,vs3,v2i,v3i)
  integer, intent(in) :: isp
  logical, intent(in) :: isperiodic
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
  real(wp), dimension(1:size(vs2,1)-4,1:size(vs2,2)-3,1:size(vs2,3)-4), intent(inout) :: v2i
  !! intent(out)
  real(wp), dimension(1:size(vs3,1)-4,1:size(vs3,2)-4,1:size(vs3,3)-3), intent(inout) :: v3i
  real(wp), dimension(-1:size(vs3,1)-2,-1:size(vs3,2)-2,-1:size(vs3,3)-2) :: param

  !> NEED TO ALSO PASS THE X2 VELOCITIES SO WE CAN COMPUTE INTERFACE VALUES
  param=vs2(:,:,:,isp)
  call halo(param,1,tag%vs2BC,isperiodic)
  !! we only need one ghost cell to compute interface velocities
  vs2(:,:,:,isp)=param

  !> PASS X3 VELOCITY BOUNDARY CONDITIONS WITH GENERIC HALOING ROUTINES
  param=vs3(:,:,:,isp)
  call halo(param,1,tag%vs3BC,isperiodic)
  !! we only need one ghost cell to compute interface velocities
  vs3(:,:,:,isp)=param
end subroutine halo_interface_vels


!> Repeat haloing operations for all species.  
subroutine halo_interface_vels_allspec(isperiodic,vs2,vs3,vs2i,vs3i,lsp)
  logical, intent(in) :: isperiodic
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
  real(wp), dimension(1:size(vs2,1)-4,1:size(vs2,2)-3,1:size(vs2,3)-4,1:size(vs2,4)), intent(inout) :: vs2i
  real(wp), dimension(1:size(vs3,1)-4,1:size(vs3,2)-4,1:size(vs3,3)-3,1:size(vs3,4)), intent(inout) :: vs3i
  integer, intent(in) :: lsp
  integer :: isp

  if (lsp>size(vs2,4)) error stop 'number of haloed species must be less than or equal to total species number'
  do isp=1,lsp
    call halo_interface_vels(isp,isperiodic,vs2,vs3,vs2i(:,:,:,isp),vs3i(:,:,:,isp))
  end do
end subroutine halo_interface_vels_allspec


!> do averaging to compute cell interface velocities; required pre-haloing in order for cells
!    to have updated boundary info from ghost cells.  Single species only, "isp"  
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


!> set values in global boundary ghost cells based on extrapolation, can be done by each worker
!    without input from the root process; does not require mpi.  This does require the inteferface
!    velocities from the x1 direction so in sense it must be preceded by an mpi call of some sort.
!    This function works for a single species "isp"
subroutine set_global_boundaries(isp,isperiodic,ns,rhovs1,vs1,vs2,vs3,rhoes,v1i)
  integer, intent(in) :: isp
  logical, intent(in) :: isperiodic
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,vs1,vs2,vs3,rhoes
  real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4), intent(in) :: v1i
  !    real(wp), parameter :: vellim=2000.0
  !    real(wp), parameter :: vellim=0.0
  real(wp) :: coeff
  integer :: ix2,ix3,lx1,lx2,lx3
  integer :: idleft,idright,idup,iddown
  real(wp), dimension(-1:size(vs3,1)-2,-1:size(vs3,2)-2,-1:size(vs3,3)-2) :: param,param2,param3,param4
  real(wp) :: tstart,tfini
  
  lx1=size(vs1,1)-4
  lx2=size(vs1,2)-4
  lx3=size(vs1,3)-4
  
  !> GHOST CELL VALUES FOR DENSITY (these need to be done last to avoid being overwritten by send/recv's!!!)
  do ix3=1,lx3
    do ix2=1,lx2
      !> logical bottom
      coeff=ns(2,ix2,ix3,isp)/ns(3,ix2,ix3,isp)
      ns(0,ix2,ix3,isp)=min(coeff*ns(1,ix2,ix3,isp),ns(1,ix2,ix3,isp))
      ns(-1,ix2,ix3,isp)=min(coeff*ns(0,ix2,ix3,isp),ns(0,ix2,ix3,isp))
  
      !> logical top
      coeff=ns(lx1-1,ix2,ix3,isp)/ns(lx1-2,ix2,ix3,isp)
      ns(lx1+1,ix2,ix3,isp)=min(coeff*ns(lx1,ix2,ix3,isp),ns(lx1,ix2,ix3,isp))
      ns(lx1+2,ix2,ix3,isp)=min(coeff*ns(lx1+1,ix2,ix3,isp),ns(lx1+1,ix2,ix3,isp))
    end do
  end do
  
  !FOR X1 MOMENTUM DENSITY
  rhovs1(0,1:lx2,1:lx3,isp)=2*v1i(1,:,:)-vs1(1,1:lx2,1:lx3,isp)
  !! initially these are velocities.  Also loose definition of 'conformable'.  Also corners never get set, but I suppose they aren't really used anyway.
  rhovs1(-1,:,:,isp)=rhovs1(0,:,:,isp)+rhovs1(0,:,:,isp)-vs1(1,:,:,isp)
  rhovs1(lx1+1,1:lx2,1:lx3,isp)=2*v1i(lx1+1,:,:)-vs1(lx1,1:lx2,1:lx3,isp)
  rhovs1(lx1+2,:,:,isp)=rhovs1(lx1+1,:,:,isp)+rhovs1(lx1+1,:,:,isp)-vs1(lx1,:,:,isp)
  
  rhovs1(-1:0,:,:,isp)=rhovs1(-1:0,:,:,isp)*ns(-1:0,:,:,isp)*ms(isp)
  !! now convert to momentum density
  rhovs1(lx1+1:lx1+2,:,:,isp)=rhovs1(lx1+1:lx1+2,:,:,isp)*ns(lx1+1:lx1+2,:,:,isp)*ms(isp)
  
  !> FOR INTERNAL ENERGY
  rhoes(0,:,:,isp)=rhoes(1,:,:,isp)
  rhoes(-1,:,:,isp)=rhoes(1,:,:,isp)
  rhoes(lx1+1,:,:,isp)=rhoes(lx1,:,:,isp)
  rhoes(lx1+2,:,:,isp)=rhoes(lx1,:,:,isp)
  
  !MZ - collect the x2 boundary conditions here - these are no longer global
  iddown=mpi_cfg%myid2-1
  idup=mpi_cfg%myid2+1
  
  !> SET THE GLOBAL X2 BOUNDARY CELLS AND ASSUME HALOING WON'T OVERWRITE.
  !> THIS DIMENSION IS ASSUMED TO NEVER BE PEREIODIC
  if (iddown==-1) then
    vs2(:,0,:,isp)=vs2(:,1,:,isp)
  
    ns(:,0,:,isp)=ns(:,1,:,isp)
    ns(:,-1,:,isp)=ns(:,1,:,isp)
    rhovs1(:,0,:,isp)=rhovs1(:,1,:,isp)
    rhovs1(:,-1,:,isp)=rhovs1(:,1,:,isp)
    rhoes(:,0,:,isp)=rhoes(:,1,:,isp)
    rhoes(:,-1,:,isp)=rhoes(:,1,:,isp)
  end if
  if (idup==mpi_cfg%lid2) then
    vs2(:,lx2+1,:,isp)=vs2(:,lx2,:,isp)
  
    ns(:,lx2+1,:,isp)=ns(:,lx2,:,isp)
    ns(:,lx2+2,:,isp)=ns(:,lx2,:,isp)
    rhovs1(:,lx2+1,:,isp)=rhovs1(:,lx2,:,isp)
    rhovs1(:,lx2+2,:,isp)=rhovs1(:,lx2,:,isp)
    rhoes(:,lx2+1,:,isp)=rhoes(:,lx2,:,isp)
    rhoes(:,lx2+2,:,isp)=rhoes(:,lx2,:,isp)
  end if
  
  !> NOW DEAL WITH ADVECTION ALONG X3; FIRST IDENTIFY MY NEIGHBORS
  idleft=mpi_cfg%myid3-1; idright=mpi_cfg%myid3+1
  
  !> SET THE GLOBAL x3 BOUNDARY CELLS AND ASSUME THAT HALOING WON'T OVERWRITE...
  if (.not. isperiodic) then
    if (idleft==-1) then
      !! left side is at global boundary, assume haloing won't overwrite
      vs3(:,:,0,isp)=vs3(:,:,1,isp)
      !! copy first cell to first ghost (vs3 not advected so only need only ghost)
  
      ns(:,:,0,isp)=ns(:,:,1,isp)
      ns(:,:,-1,isp)=ns(:,:,1,isp)
      rhovs1(:,:,0,isp)=rhovs1(:,:,1,isp)
      rhovs1(:,:,-1,isp)=rhovs1(:,:,1,isp)
      rhoes(:,:,0,isp)=rhoes(:,:,1,isp)
      rhoes(:,:,-1,isp)=rhoes(:,:,1,isp)
    end if
    if (idright==mpi_cfg%lid3) then    !my right boundary is the global boundary, assume haloing won't overwrite
      vs3(:,:,lx3+1,isp)=vs3(:,:,lx3,isp)    !copy last cell to first ghost (all that's needed since vs3 not advected)
  
      ns(:,:,lx3+1,isp)=ns(:,:,lx3,isp)
      ns(:,:,lx3+2,isp)=ns(:,:,lx3,isp)
      rhovs1(:,:,lx3+1,isp)=rhovs1(:,:,lx3,isp)
      rhovs1(:,:,lx3+2,isp)=rhovs1(:,:,lx3,isp)
      rhoes(:,:,lx3+1,isp)=rhoes(:,:,lx3,isp)
      rhoes(:,:,lx3+2,isp)=rhoes(:,:,lx3,isp)
    end if
  end if
end subroutine set_global_boundaries


!> set global boundaries for all species
subroutine set_global_boundaries_allspec(isperiodic,ns,rhovs1,vs1,vs2,vs3,rhoes,vs1i,lsp)
  logical, intent(in) :: isperiodic
  real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: ns,rhovs1,vs1,vs2,vs3,rhoes
  real(wp), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4,1:size(vs1,4)), intent(in) :: vs1i
  integer, intent(in) :: lsp
  integer :: isp

  if (lsp>size(vs1,4)) error stop 'number of global boundaries must be less than or equal to total species number'
  do isp=1,lsp
    call set_global_boundaries(isp,isperiodic,ns,rhovs1,vs1,vs2,vs3,rhoes,vs1i(:,:,:,isp))
  end do
end subroutine set_global_boundaries_allspec


!> 3D advection subroutine
function advec3D_MC_mpi_curv_23(f,v1i,v2i,v3i,dt,x,frank,tagf)
!------------------------------------------------------------
!-------ADVECT A VARIABLE IN 3D FOR AN MPI SIMULATION
!------------------------------------------------------------
!-------It is critical that the mpi'd dimension be advected
!-------first to avoid having to repass ghost cells between
!-------workers after 1 and 2 dimension are advected.
!-------
!-------NOTE: also that the ghost cells should really
!-------be updated after each sweep.  Ie the x2 boundary regions
!-------should be updated after x3 sweep and the x1 boundary
!-------conditions should be updated after x3,x2 sweeps.  I'm
!-------really to lazy to deal with this now...
  real(wp), dimension(-1:,-1:,-1:), intent(in) :: f    !f includes ghost cells
  real(wp), dimension(:,:,:), intent(in) :: v1i
  real(wp), dimension(:,:,:), intent(in) :: v2i
  real(wp), dimension(:,:,:), intent(in) :: v3i
  real(wp), intent(in) :: dt
  class(curvmesh), intent(in) :: x
  integer, intent(in) :: frank    !f's rank so that we know which metric coeffs to use.
  integer, intent(in) :: tagf
  real(wp), dimension(-1:size(f,1)-2,-1:size(f,2)-2,-1:size(f,3)-2) :: advec3D_MC_mpi_curv_23
  
  !We are assuming here that the data have been pre-haloed before this function is called.  This must be the case
  !to avoid tearing artifacts...  In the future this should be removed in favor of explicit haloing withint
  !this procedure, which is already required for x2 message passsing.
  
  !we must recopying the boundary conditions after any halo operation (which assumes periodic)
  advec3D_MC_mpi_curv_23=f
  call halo(advec3D_MC_mpi_curv_23,2,tagf,x%flagper)
  
  !x3-sweep
  call sweep3(advec3D_MC_mpi_curv_23,v3i,dt,x,frank)
  
  !x1-sweep
  call sweep1(advec3D_MC_mpi_curv_23,v1i,dt,x)
  
  !at this point if we've divided in two dimensions with mpi it is necessary to halo again before
  !the final sweep to avoid tearing artifacts...
  call halo(advec3D_MC_mpi_curv_23,2,tagf,x%flagper)
  
  !x2-sweep, if necessary
  call sweep2(advec3D_MC_mpi_curv_23,v2i,dt,x,frank)
end function advec3D_MC_mpi_curv_23


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

end module advec_mpi
