submodule (sources) sources_mpi

use mpimod, only: myid, tagvs1bc, tagvs2bc, tagvs3bc, lid, halo, myid2,myid3,lid2,lid3
implicit none

contains

module procedure RK2_prep_mpi
!subroutine RK2_prep_mpi(isp,isperiodic,vs1,vs2,vs3)
!! PASS BOUNDARY CELLS FOR COMPUTING COMPRESSION.
!! DONE ON A PER-SPECIES BASIS.
!! ION PARAMETER ARGUMENTS SHOULD INCLUDE GHOST CELLS.
!! DO WE NEED TO PASS V1,2 VARIABLES FOR DIV?

real(wp), dimension(-1:size(vs1,1)-2,-1:size(vs1,2)-2,-1:size(vs1,3)-2) :: param
integer :: lx1,lx2,lx3
integer :: idleft,idright,idup,iddown


lx1=size(vs1,1)-4
lx2=size(vs1,2)-4
lx3=size(vs1,3)-4


!ZOH EXTRAPOLATION OF V1,2 VARIABLES
vs1(0,:,:,isp)=vs1(1,:,:,isp)
vs1(lx1+1,:,:,isp)=vs1(lx1,:,:,isp)


!IDENTIFY MY NEIGHBORS in x2 and x3
idleft=myid3-1; idright=myid3+1
iddown=myid2-1; idup=myid2+1

!-- Now halo the interior parts (must happen for every worker since even a worker with a
!-- global boundary will still have one interior boundary to be haloed.
!BY DEFAULT THE GLOBAL BOUNDARIES ARE ASSUMED TO BE PERIOIDIC
param=vs1(:,:,:,isp)
call halo(param,1,tagvs1BC,isperiodic)
vs1(:,:,:,isp)=param
param=vs2(:,:,:,isp)
call halo(param,1,tagvs2BC,isperiodic)
vs2(:,:,:,isp)=param
param=vs3(:,:,:,isp)
call halo(param,1,tagvs3BC,isperiodic)
vs3(:,:,:,isp)=param


!ZERO ORDER HOLD EXTRAPOLATION OF BOUNDARIES (UNLESS PERIODIC)
if(iddown==-1) then
  vs1(:,0,:,isp)=vs1(:,1,:,isp)
  vs2(:,0,:,isp)=vs2(:,1,:,isp)
  vs3(:,0,:,isp)=vs3(:,1,:,isp)
end if
if(idup==lid2) then
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
  if (idright==lid3) then    !right x3 boundary
    vs1(:,:,lx3+1,isp)=vs1(:,:,lx3,isp)
    vs2(:,:,lx3+1,isp)=vs2(:,:,lx3,isp)
    vs3(:,:,lx3+1,isp)=vs3(:,:,lx3,isp)
  end if
end if

end procedure RK2_prep_mpi

end submodule sources_mpi