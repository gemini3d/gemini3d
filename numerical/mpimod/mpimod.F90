module mpimod
!! NOTES:
!! * Need to consider overloading routines as send_ghost and send_noghost so that
!!   it is more clear what the structure of the input arrays should be.

use phys_consts, only : lsp, wp
!! code needs to know how many species are being used.

use mpi, only: mpi_init, mpi_comm_world, &
               mpi_integer,mpi_sum, &
               mpi_status_size, mpi_status_ignore, &
#if REALBITS==32
mpi_realprec=>mpi_real
#elif REALBITS==64
mpi_realprec=>mpi_double_precision
#endif

implicit none

private :: lsp, wp
private :: mpi_init, mpi_status_size, mpi_status_ignore, mpi_integer, mpi_sum, mpi_comm_world

!> NOW A LIST OF TAGS SO THESE DO NOT NEED TO BE EMBEDDED IN EACH SUBROUTINE
integer, parameter :: tagns=2,tagvs1=3,tagTs=4
!! root/workers input routines.  also output routines for root/worker

integer, parameter :: tagJ1=6,tagJ2=7,tagJ3=8,tagv2=9,tagv3=10
!! output root/worker routines, main program, potential_comm
integer, parameter :: tagvs3BC=100,tagnsBC=101,tagrhovs1BC=102,tagrhoesBC=103
!! used in the advection boundary conditions
integer, parameter :: tagvs1BC=1000,tagvs2BC=1001!,tagvs3BC=1002    !used in the compression solution

!> THESE MESSAGES ARE USED IN ELECTRODYNAMICS MODULE
integer, parameter :: tagE1=11,tagE2=12,tagE3=13
integer, parameter :: tagsigP=16,tagsigH=17,tagsig0=18,tagincap=19,tagv2pol=20,tagv3pol=21
integer, parameter :: tagDE2Dt=22,tagDE3Dt=23,tagflagdirich=24    !unused in present version
integer, parameter :: tagvn2=25,tagvn3=26,tagB1=27                !for passing/gathering full-grid winds

!> IN THE MAIN PROGRAM
integer, parameter :: tagx3=1,tagdt=5
integer, parameter :: tagx1=27,tagx2=28

!> IN THE GRID MODULE
integer, parameter :: tagh1=29,tagh2=30,tagh3=31
integer, parameter :: tagglat=32,tagglon=33,tagalt=34
integer, parameter :: taglx1=35,taglx2=36,taglx3=37,taglx3all=38
integer, parameter :: tagBmag=39,taginc=40,tagnull=41
integer, parameter :: tageunit1=42,tageunit2=43,tageunit3=44,tager=45,tagetheta=46,tagephi=47
integer, parameter :: tagr=56,tagtheta=57,tagphi=58

!> IN THE NEUTRAL MODULE
integer, parameter :: taglrho=48,taglz=49
integer, parameter :: tagdnO=50,tagdnN2=51,tagdnO2=52,tagdTn=53,tagdvnrho=54,tagdvnz=55
integer, parameter :: tagly=69

!> FOR DEALING WITH PRECIPITATION BOUNDARY CONDITIONS MODULE
integer, parameter :: tagllat=59,tagllon=60,tagmlat=61,tagmlon=62,tagQp=63,tagE0p=64

!> FOR DEALING WITH THE ELECTRIC FIELD BOUNDARY CONDITIONS
integer, parameter :: tagE0xp=65,tagE0yp=66,tagE0xi=67,tagE0yi=68

!> FOR DISTRIBUTING PART OF THE ELECTRODYNAMICS CALCULATIONS
integer, parameter :: tagsrc=69,tagSigPint2=70,tagSigPint3=71,tagSigHint=72,tagincapint=73,tagv2electro=74,tagv3electro=75
integer, parameter :: tagE01=76,tagE02=77,tagE03=78,tagVminx1=79,tagVmaxx1=80

!> THESE ARE USED IN MAGCALC.F90 PROGRAM
integer, parameter :: tagBr=81,tagBtheta=82, tagBphi=83
integer, parameter :: tagdV=84,tagJx=85,tagJy=86,tagRx=87,tagRy=88,tagRz=89,tagRcubed=90,tagJz=91

!> FOR COMMUNICATING IF THE GRID DIMENSIONS HAVE BEEN SWAPPED
integer, parameter :: tagswap=92

!> AURORAL TAG(S)
integer, parameter :: tagAur=95

!> VARIABLES REUSED BY ALL WORKERS AND USING MODULES
integer, protected :: myid,lid
!! no external procedure should mess with these (but they need to be able to read them)

integer :: ierr
!! using procedures need to be able to overwrite this to prevent seg. faults (or something)


!> VARIABLES RELATED TO PROCESS GRID (IF USED)
integer, protected :: lid2,lid3,myid2,myid3


!> THESE INTERFACES OVERLOAD THE MPI GATHER,BROADCAST SUBROUTINES FOR ARRAYS OF DIFFERENT RANKS.
interface gather_recv
  procedure gather_recv2D, gather_recv3D, gather_recv4D
end interface gather_recv

interface gather_send
  procedure gather_send2D, gather_send3D, gather_send4D
end interface gather_send

interface ! mpisend

module subroutine gather_send2D(paramtrim,tag)
real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
end subroutine gather_send2D

module subroutine gather_send3D(paramtrim,tag)
real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
end subroutine gather_send3D

module subroutine gather_send4D(param,tag)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
end subroutine gather_send4D

module subroutine bcast_send1D(paramall,tag,param)
real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param
end subroutine bcast_send1D

module subroutine bcast_send2D(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrim
end subroutine bcast_send2D

module subroutine bcast_send3D(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim
end subroutine bcast_send3D

module subroutine bcast_send4D(paramall,tag,param)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
end subroutine bcast_send4D

module subroutine bcast_send3D_x3i(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim
end subroutine bcast_send3D_x3i

module subroutine bcast_send3D_ghost(paramall,tag,param)
real(wp), dimension(-1:,-1:,-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
end subroutine bcast_send3D_ghost

end interface

interface ! mpirecv

module subroutine gather_recv2D(paramtrim,tag,paramtrimall)
real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrimall
end subroutine gather_recv2D

module subroutine gather_recv3D(paramtrim,tag,paramtrimall)
real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrimall
end subroutine gather_recv3D

module subroutine gather_recv4D(param,tag,paramall)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: paramall
end subroutine gather_recv4D

module subroutine bcast_recv1D(param,tag)
real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv1D

module subroutine bcast_recv2D(paramtrim,tag)
real(wp), dimension(:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv2D

module subroutine bcast_recv3D(paramtrim,tag)
real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv3D

module subroutine bcast_recv4D(param,tag)
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv4D

module subroutine bcast_recv3D_x3i(paramtrim,tag)
real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv3D_x3i

module subroutine bcast_recv3D_ghost(param,tag)
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv3D_ghost

end interface

interface ! mpihalo

module subroutine halo(param,lhalo,tag)
real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (1 or 2 only)
integer, intent(in) :: tag
end subroutine halo

module subroutine halo_end(param,paramend,tag)
real(wp), dimension(:,:,:), intent(inout) :: param
real(wp), dimension(:,:), intent(out) :: paramend
integer, intent(in) :: tag
end subroutine halo_end

end interface

interface bcast_send
  procedure bcast_send1D, bcast_send2D, bcast_send3D, bcast_send4D
end interface bcast_send

interface bcast_recv
  procedure bcast_recv1D, bcast_recv2D, bcast_recv3D, bcast_recv4D
end interface bcast_recv


contains


subroutine mpisetup()
!! INITIALIZES MODULE MPI VARIABLES FOR A WORKER.
!! THIS CURRENTLY IS UNUSED AS IT HAS NOT BEEN
!! FULLY IMPLEMENTED IN THIS VERSINO OF THE CODE.

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,lid,ierr)

!> INITIALIZE, ONLY PARALLELIZING IN X3, GRIDDING FUNCTION MAY CHANGE THIS, IF CALLED.
  lid2=1
  lid3=lid

end subroutine mpisetup


subroutine mpigrid(lx2all,lx3all)
!! THIS SUBROUTINE DEFINES A PROCESS GRID, IF REQUIRED
!! IT IS CURRENTLY NOT USED BUT KEPT HERE FOR FUTURE DEVELOPMENT

integer, intent(in) :: lx2all,lx3all

lid2=1
lid3=lid
do while( ((lid3/2)*2==lid3) .and. (lid3-lid2>lid3 .or. lid3-lid2>lid2) .and. &
         lx3all/(lid3/2)*(lid3/2)==lx3all .and. lx2all/(lid2*2)*(lid2*2)==lx2all)
!! must insure that lx3 is divisible by lid3 and lx2 by lid2

  lid3=lid3/2
  lid2=lid2*2
end do


!THIS PROCESS' LOCATION ON THE GRID
myid3=myid/lid2
myid2=myid-myid3*lid2

print *, 'Proposed process grid is x2 by x3 size (in number of processes):  ',lid2,' by ',lid3
print *, 'Process:  ',myid,' is at location:  ',myid2,myid3,' on the process grid'

end subroutine mpigrid


subroutine mpibreakdown()
!! SHUTS DOWN MPI

call mpi_finalize(ierr)

end subroutine mpibreakdown


end module mpimod
