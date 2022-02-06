module advec_mpi

use phys_consts, only: lsp,ms, wp
use grid, only : gridflag
use meshobj, only: curvmesh
  !! do not import grid sizes in case we want do subgrid advection...
use mpimod, only: mpi_cfg, halo, tag=>gemini_mpi

implicit none (type, external)
private
public :: halo_interface_vels_allspec,set_global_boundaries_allspec
contains
  !> Perform haloing needed to ghost-fill so cell interface vels (single species) can be computed across the grid
  subroutine halo_interface_vels(isp,isperiodic,vs2,vs3)
    integer, intent(in) :: isp
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
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
  subroutine halo_interface_vels_allspec(isperiodic,vs2,vs3,lsp)
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs2,vs3
    integer, intent(in) :: lsp
    integer :: isp
  
    if (lsp>size(vs2,4)) error stop 'number of haloed species must be less than or equal to total species number'
    do isp=1,lsp
      call halo_interface_vels(isp,isperiodic,vs2,vs3)
    end do
  end subroutine halo_interface_vels_allspec
  
  
  !> set values in global boundary ghost cells based on extrapolation, can be done by each worker
  !    without input from the root process; does not require mpi.  This does require the inteferface
  !    velocities from the x1 direction so in sense it must be preceded by an mpi call of some sort.
  !    This function works for a single species "isp".  It does need information about the mpi image
  !    organization and as such is considered part of the advec_mpi module.
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
end module advec_mpi
