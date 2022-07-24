!> This module contains mpi-related subroutines needed for solving source-loss numerical problems in GEMINI
module sources_mpi

use phys_consts, only: wp
use mpimod, only: mpi_cfg, tag=>gemini_mpi, halo

implicit none (type, external)
private
public :: RK2_prep_mpi_allspec, RK2_global_boundary_allspec

contains
  !> This haloes a single ghost cell for just the three components of velocity so a divergence
  !    can be calculated.  
  subroutine RK2_prep_mpi(isp,isperiodic,vs1,vs2,vs3)
    integer, intent(in) :: isp
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    !subroutine RK2_prep_mpi(isp,isperiodic,vs1,vs2,vs3)
    !! PASS BOUNDARY CELLS FOR COMPUTING COMPRESSION.
    !! DONE ON A PER-SPECIES BASIS.
    !! ION PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS.
    !! DO WE NEED TO PASS V1,2 VARIABLES FOR DIV?
    real(wp), dimension(-1:size(vs1,1)-2,-1:size(vs1,2)-2,-1:size(vs1,3)-2) :: param
            
    !-- Now halo the interior parts (must happen for every worker since even a worker with a
    !-- global boundary will still have one interior boundary to be haloed.
    !BY DEFAULT THE GLOBAL BOUNDARIES ARE ASSUMED TO BE PERIOIDIC
    param=vs1(:,:,:,isp)
    call halo(param,1,tag%vs1BC,isperiodic)
    vs1(:,:,:,isp)=param
    param=vs2(:,:,:,isp)
    call halo(param,1,tag%vs2BC,isperiodic)
    vs2(:,:,:,isp)=param
    param=vs3(:,:,:,isp)
    call halo(param,1,tag%vs3BC,isperiodic)
    vs3(:,:,:,isp)=param
  end subroutine RK2_prep_mpi


  !> halo all species velocities in order to be ready for compression substep
  subroutine RK2_prep_mpi_allspec(vs1,vs2,vs3,isperiodic)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    logical, intent(in) :: isperiodic
    integer isp,lsp
  
    lsp=size(vs1,4)
    do isp=1,lsp
      call RK2_prep_mpi(isp,isperiodic,vs1,vs2,vs3)    !role-agnostic mpi, all-to-neighbor
    end do
    !call RK2_global_boundary_allspec(vs1,vs2,vs3,isperiodic)     ! separate call now
  end subroutine RK2_prep_mpi_allspec


  !> global boundaries for all species
  subroutine RK2_global_boundary_allspec(vs1,vs2,vs3,isperiodic)
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    logical, intent(in) :: isperiodic
    integer isp,lsp

    lsp=size(vs1,4)
    do isp=1,lsp
      call RK2_global_boundary(isp,isperiodic,vs1,vs2,vs3)
    end do
  end subroutine RK2_global_boundary_allspec


  !> correct/extrapolate global boundaries for compression substep
  subroutine RK2_global_boundary(isp,isperiodic,vs1,vs2,vs3)
    integer, intent(in) :: isp
    logical, intent(in) :: isperiodic
    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: vs1,vs2,vs3
    integer :: lx1,lx2,lx3
    integer :: idright,idleft,idup,iddown

    ! convenience sizes
    lx1=size(vs1,1)-4
    lx2=size(vs1,2)-4
    lx3=size(vs1,3)-4

    !IDENTIFY MY NEIGHBORS in x2 and x3
    idleft=mpi_cfg%myid3-1; idright=mpi_cfg%myid3+1
    iddown=mpi_cfg%myid2-1; idup=mpi_cfg%myid2+1

    !ZOH EXTRAPOLATION OF V1 VARIABLES in x1
    vs1(0,:,:,isp)=vs1(1,:,:,isp)
    vs1(lx1+1,:,:,isp)=vs1(lx1,:,:,isp)

    !ZERO ORDER HOLD EXTRAPOLATION OF BOUNDARIES (UNLESS PERIODIC)
    if(iddown==-1) then
      vs1(:,0,:,isp)=vs1(:,1,:,isp)
      vs2(:,0,:,isp)=vs2(:,1,:,isp)
      vs3(:,0,:,isp)=vs3(:,1,:,isp)
    end if
    if(idup==mpi_cfg%lid2) then
      vs1(:,lx2+1,:,isp)=vs1(:,lx2,:,isp)
      vs2(:,lx2+1,:,isp)=vs2(:,lx2,:,isp)
      vs3(:,lx2+1,:,isp)=vs3(:,lx2,:,isp)
    end if
    if (.not. isperiodic) then
      if (idleft==-1) then    !left x3 boundary
        vs1(:,:,0,isp)=vs1(:,:,1,isp)
        vs2(:,:,0,isp)=vs2(:,:,1,isp)
        vs3(:,:,0,isp)=vs3(:,:,1,isp)
      end if
      if (idright==mpi_cfg%lid3) then    !right x3 boundary
        vs1(:,:,lx3+1,isp)=vs1(:,:,lx3,isp)
        vs2(:,:,lx3+1,isp)=vs2(:,:,lx3,isp)
        vs3(:,:,lx3+1,isp)=vs3(:,:,lx3,isp)
      end if
    end if 
  end subroutine RK2_global_boundary
end module sources_mpi
