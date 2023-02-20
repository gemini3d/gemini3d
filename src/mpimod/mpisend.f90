submodule (mpimod) mpisend

use mpi_f08, only : mpi_send

implicit none (type, external)

contains


!subroutine gather_send2D_23(paramtrim,tag)
!real(wp), dimension(:,:), intent(in) :: paramtrim
!integer, intent(in) :: tag
module procedure gather_send2D_23
!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
!-------ANY GHOST CELLS!!!!
!-------THIS ROUTINE WORKS ON A PROCESS GRID
!------------------------------------------------------------


integer :: lx2,lx3


lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx3=size(paramtrim,2)

call mpi_send(paramtrim, lx2*lx3, mpi_realprec, 0, tag, MPI_COMM_WORLD)

end procedure gather_send2D_23


!subroutine gather_send3D_23(paramtrim,tag)
!real(wp), dimension(:,:,:), intent(in) :: paramtrim
!integer, intent(in) :: tag
module procedure gather_send3D_23
!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------ANY GHOST CELLS!!!!
!-------THIS VERSION WORKS ON A PROCESS GRID
!------------------------------------------------------------


integer :: lx1,lx2,lx3


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)

call mpi_send(paramtrim,lx1*lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD)

end procedure gather_send3D_23


!subroutine gather_send4D_23(param,tag)
!real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
!integer, intent(in) :: tag
module procedure gather_send4D_23
!------------------------------------------------------------
!-------SENDS 4D DATA ON A 2D PROCESS GRID TO ROOT.
!------------------------------------------------------------


integer :: lx1,lx2,lx3,isp
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4

do isp=1,lsp
  paramtmp=param(-1:lx1+2,1:lx2,1:lx3,isp)
  call mpi_send(paramtmp,(lx1+4)*lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD)
end do

end procedure gather_send4D_23


!> send our data to root, ***including*** ghost cells, note that the ghost cells
!    need to be populated correctly or this will stripe the data.
module procedure gather_send3D_ghost_23
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!!ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH INCLUDE GHOST CELLS

  integer :: lx1,lx2,lx3

  !> note here that param has ghost cells
  lx1=size(param,1)-4
  lx2=size(param,2)-4
  lx3=size(param,3)-4

  !> workers send their slab of data to root
  call mpi_send(param,(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,0,tag,MPI_COMM_WORLD)
end procedure gather_send3D_ghost_23


module procedure gather_send3D_x2i_23

  integer :: lx1,lx2,lx3

  !> note here that param has ghost cells
  lx1=size(param,1)
  lx2=size(param,2)-1     ! this dim has extra interface val.
  lx3=size(param,3)

  !> workers send their slab of data to root
  call mpi_send(param,(lx1)*(lx2+1)*(lx3), &
               mpi_realprec,0,tag,MPI_COMM_WORLD)
end procedure gather_send3D_x2i_23


module procedure gather_send3D_x3i_23

  integer :: lx1,lx2,lx3

  !> note here that param has ghost cells
  lx1=size(param,1)
  lx2=size(param,2)
  lx3=size(param,3)-1

  !> workers send their slab of data to root
  call mpi_send(param,(lx1)*(lx2)*(lx3+1), &
               mpi_realprec,0,tag,MPI_COMM_WORLD)
end procedure gather_send3D_x3i_23


!subroutine bcast_send1D_23_2(paramall,tag,param)
!real(wp), dimension(-1:), intent(in) :: paramall
!integer, intent(in) :: tag
!real(wp), dimension(-1:), intent(out) :: param
module procedure bcast_send1D_23_2
!------------------------------------------------------------
!-------BROADCASTS MPI DIMENSION VARIABLES TO WORKERS.  NOTE THAT
!-------WE'VE ELECTED TO NOT USE THE GENERAL BROADCAST ROUTINES FOR
!-------SINCE THESE OPERATIONS REQUIRE A LOT OF SPECIAL CASING FOR
!-------THE SIZES OF THE VARIABLES TO BE SENT
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 1D ARRAYS THAT ARE REPRESENTING
!-------THE X2-DIRECTION
!------------------------------------------------------------


integer :: lx,lxall     !local sizes
integer :: iid,islstart,islfin
integer, dimension(2) :: indsgrid

lxall=size(paramall,1)-4
lx=size(param,1)-4


do iid=1,mpi_cfg%lid-1
  indsgrid=ID2grid(iid, mpi_cfg%lid2)

  islstart=indsgrid(1)*lx+1
  islfin=islstart+lx-1

  call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
end do
param=paramall(-1:lx+2)

end procedure bcast_send1D_23_2


!subroutine bcast_send1D_23_3(paramall,tag,param)
!real(wp), dimension(-1:), intent(in) :: paramall
!integer, intent(in) :: tag
!real(wp), dimension(-1:), intent(out) :: param
module procedure bcast_send1D_23_3
!------------------------------------------------------------
!-------BROADCASTS MPI DIMENSION VARIABLES TO WORKERS.  NOTE THAT
!-------WE'VE ELECTED TO NOT USE THE GENERAL BROADCAST ROUTINES FOR
!-------SINCE THESE OPERATIONS REQUIRE A LOT OF SPECIAL CASING FOR
!-------THE SIZES OF THE VARIABLES TO BE SENT
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 1D ARRAYS THAT ARE REPRESENTING
!-------THE X3-DIRECTION
!------------------------------------------------------------


integer :: lx,lxall     !local sizes
integer :: iid,islstart,islfin
integer, dimension(2) :: indsgrid

lxall=size(paramall,1)-4
lx=size(param,1)-4


do iid=1,mpi_cfg%lid-1
  indsgrid=ID2grid(iid, mpi_cfg%lid2)    !compute my location on the process grid

  islstart=indsgrid(2)*lx+1    !send piece of grid that corresponds to my x3 position
  islfin=islstart+lx-1

  call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
end do
param=paramall(-1:lx+2)

end procedure bcast_send1D_23_3


!subroutine bcast_send2D_23(paramtrimall,tag,paramtrim)
!real(wp), dimension(:,:), intent(in) :: paramtrimall
!integer, intent(in) :: tag
!real(wp), dimension(:,:), intent(out) :: paramtrim
module procedure bcast_send2D_23
!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------


integer :: lx2,lx3
integer :: iid
integer, dimension(4) :: inds

real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2)) :: paramtmp


lx2=size(paramtrim,1)    !assume this is an array which has been 'flattened' along the 1-dimension
lx3=size(paramtrim,2)


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,mpi_cfg%lid-1
!  islstart=iid*lx3+1
!  islfin=islstart+lx3-1
  inds=slabinds(iid,lx2,lx3)

  paramtmp=paramtrimall(inds(1):inds(2),inds(3):inds(4))
  call mpi_send(paramtmp,lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim(1:lx2,1:lx3)=paramtrimall(1:lx2,1:lx3)

end procedure bcast_send2D_23


!subroutine bcast_send3D_23(paramtrimall,tag,paramtrim)
!real(wp), dimension(:,:,:), intent(in) :: paramtrimall
!integer, intent(in) :: tag
!real(wp), dimension(:,:,:), intent(out) :: paramtrim
module procedure bcast_send3D_23
!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS!
!-------
!-------ALSO NOTE THAT IF THE ARRAY SIZE (DIM 3)  DOES NOT CORRESPOND
!-------TO THE SIZE OF THE SYSTEM IN THE X3-DIRECTION, THEN
!-------THE SLAB CALCULATIONS FOR WORKERS WILL BE OFF.
!------------------------------------------------------------


integer :: lx1,lx2,lx3
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,mpi_cfg%lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(1:lx1,inds(1):inds(2),inds(3):inds(4))
  call mpi_send(paramtmp,lx1*lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
end do


!> ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(1:lx1,1:lx2,1:lx3)

end procedure bcast_send3D_23


!subroutine bcast_send3D_x3i_23(paramtrimall,tag,paramtrim)
!real(wp), dimension(:,:,:), intent(in) :: paramtrimall
!integer, intent(in) :: tag
!real(wp), dimension(:,:,:), intent(out) :: paramtrim
module procedure bcast_send3D_x3i_23
!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS, BUT ARE X3 INTERFACE QUANITTIES AND HENCE
!-------LARGER THAN  LX3
!------------------------------------------------------------


integer :: lx1,lx2,lx3
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp
!! has size lx3+1 due to input having that size

lx1=size(paramtrim,1)      !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)-1    !note that we are interpreting input as an x3i quantity meaning that it has size lx3+1


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,mpi_cfg%lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(:,inds(1):inds(2),inds(3):inds(4)+1)
  !! +1 since this is an x3 interface quantity
  call mpi_send(paramtmp,lx1*lx2*(lx3+1), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
  !! note the +1 since these are interface quantities (and need to overlap b/t workers)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,1:lx2,1:lx3+1)

end procedure bcast_send3D_x3i_23


!subroutine bcast_send3D_x2i_23(paramtrimall,tag,paramtrim)
!real(wp), dimension(:,:,:), intent(in) :: paramtrimall
!integer, intent(in) :: tag
!real(wp), dimension(:,:,:), intent(out) :: paramtrim
module procedure bcast_send3D_x2i_23
!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS, BUT ARE X3 INTERFACE QUANITTIES AND HENCE
!-------LARGER THAN  LX2
!------------------------------------------------------------


integer :: lx1,lx2,lx3
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp
!! has size lx3+1 due to input having that size

lx1=size(paramtrim,1)      !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)-1    !for this routine the interface is in the x2-direciton
lx3=size(paramtrim,3)


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,mpi_cfg%lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(:,inds(1):inds(2)+1,inds(3):inds(4))
  !! +1 since this is an x3 interface quantity
  call mpi_send(paramtmp,lx1*(lx2+1)*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
  !! note the +1 since these are interface quantities (and need to overlap b/t workers)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,1:lx2+1,1:lx3)

end procedure bcast_send3D_x2i_23


!subroutine bcast_send3D_ghost_23(paramall,tag,param)
!real(wp), dimension(-1:,-1:,-1:), intent(in) :: paramall
!integer, intent(in) :: tag
!real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
module procedure bcast_send3D_ghost_23
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!!ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH INCLUDE GHOST CELLS


  integer :: lx1,lx2,lx3
  integer :: iid
  integer, dimension(4) :: inds
  real(wp), dimension(-1:size(param,1)-2,-1:size(param,2)-2,-1:size(param,3)-2) :: paramtmp

  !> note here that param has ghost cells
  lx1=size(param,1)-4
  lx2=size(param,2)-4
  lx3=size(param,3)-4

  !> ROOT BROADCASTS IC DATA TO WORKERS
  do iid=1,mpi_cfg%lid-1
    inds=slabinds(iid,lx2,lx3)
    paramtmp=paramall(:,inds(1)-2:inds(2)+2,inds(3)-2:inds(4)+2)
    call mpi_send(paramtmp,(lx1+4)*(lx2+4)*(lx3+4), &
                 mpi_realprec,iid,tag,MPI_COMM_WORLD)
  end do

  !> ROOT TAKES A SLAB OF DATA
  param=paramall(:,-1:lx2+2,-1:lx3+2)
end procedure bcast_send3D_ghost_23


!subroutine bcast_send4D_23(paramall,tag,param)
!real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
!integer, intent(in) :: tag
!real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
module procedure bcast_send4D_23
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!! ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!! GHOST CELLS!


integer :: lx1,lx2,lx3,isp
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> ROOT BROADCASTS IC DATA TO WORKERS
do isp=1,lsp
  param(:,1:lx2,1:lx3,isp)=paramall(:,1:lx2,1:lx3,isp)
  !! roots part of the data

  do iid=1,mpi_cfg%lid-1
    inds=slabinds(iid,lx2,lx3)
    paramtmp=paramall(-1:lx1+2,inds(1):inds(2),inds(3):inds(4),isp)
    call mpi_send(paramtmp,(lx1+4)*lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD)
  end do
end do

end procedure bcast_send4D_23


end submodule mpisend
