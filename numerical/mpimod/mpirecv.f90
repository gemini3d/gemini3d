submodule (mpimod) mpirecv

implicit none

contains

module procedure gather_recv2D
!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!!
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!!
!! THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE ANY GHOST CELLS!!!!

integer :: lx1,lx2,lx3,isp,lsp
integer :: iid,islstart,islfin


lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx3=size(paramtrim,2)


!PATCH DATA TOGETHER FOR OUTPUT STARTING WITH ROOT'S SLAB
paramtrimall(:,1:lx3)=paramtrim   !copy root's data into full-grid array

do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_recv(paramtrimall(:,islstart:islfin),lx2*lx3, &
                mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end do

end procedure gather_recv2D


module procedure gather_recv3D
!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!!
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! ANY GHOST CELLS!!!!

integer :: lx1,lx2,lx3,isp,lsp
integer :: iid,islstart,islfin


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!PATCH DATA TOGETHER FOR OUTPUT STARTING WITH ROOT'S SLAB
paramtrimall(:,:,1:lx3)=paramtrim   !copy root's data into full-grid array

do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_recv(paramtrimall(:,:,islstart:islfin),lx1*lx2*lx3, &
                mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end do

end procedure gather_recv3D


module procedure gather_recv4D

!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

integer :: lx1,lx2,lx3,isp
integer :: iid,islstart,islfin

lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!Originally the outer loop was over worker number, which cycles the 3rd dimension
!slower than 4th (backward from what would be most efficient memory access pattern)
!Since the gathering operation is root-limited probably, I'm guessing it's better
!to give root an efficient memory access pattern here, but I haven't tested this
!theory.
do isp=1,lsp
  paramall(:,:,1:lx3,isp)=param(:,:,1:lx3,isp)

  do iid=1,lid-1
    islstart=iid*lx3+1
    islfin=islstart+lx3-1

    call mpi_recv(paramall(:,:,islstart:islfin,isp),(lx1+4)*(lx2+4)*lx3, &
                  mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  end do
end do

end procedure gather_recv4D


module procedure bcast_recv1D
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 1D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

integer :: lx
integer :: iid

lx=size(param,1)-4

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx+4), &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end procedure bcast_recv1D


module procedure bcast_recv2D
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

integer :: lx2,lx3
integer :: iid


lx2=size(paramtrim,1)
lx3=size(paramtrim,2)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx2*lx3, &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end procedure bcast_recv2D


module procedure bcast_recv3D
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

integer :: lx1,lx2,lx3
integer :: iid

!> note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*lx2*lx3, &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end procedure bcast_recv3D


module procedure bcast_recv3D_x3i
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

integer :: lx1,lx2,lx3
integer :: iid

!>note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)-1
  !! `lx3` is an x3i quantity

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*lx2*(lx3+1), &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end procedure bcast_recv3D_x3i


module procedure bcast_recv3D_ghost
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

integer :: lx1,lx2,lx3
integer :: iid

!> note here that param has ghost cells
lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end procedure bcast_recv3D_ghost


module procedure bcast_recv4D
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

integer :: lx1,lx2,lx3,isp

lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!WORKERS RECEIVE THE IC DATA FROM ROOT
do isp=1,lsp
  call mpi_recv(param(:,:,:,isp),(lx1+4)*(lx2+4)*(lx3+4), &
                 mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end do

end procedure bcast_recv4D


end submodule mpirecv
