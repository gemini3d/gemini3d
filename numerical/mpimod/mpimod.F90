module mpimod
!! NOTES:
!! * Need to consider overloading routines as send_ghost and send_noghost so that
!!   it is more clear what the structure of the input arrays should be.
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only : lsp, wp
!! code needs to know how many species are being used.

use mpi, only: mpi_init, mpi_comm_world, &
               mpi_integer,mpi_sum, &
               mpi_status_size, mpi_status_ignore, MPI_PROC_NULL, &

#if REALBITS==32
mpi_realprec=>mpi_real
#else
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

!> FOR SENDING THE FULL X2 GRID SIZE
integer, parameter :: taglx2all=93
integer, parameter :: tagx2all=94
integer, parameter :: tagx3all=95

!> AURORAL TAG(S)
integer, parameter :: tagAur=96, tagZxden=200

!!> GENERIC PARAMETER (USED BY ADVECTION CODE - HOPEFULLY DOESN'T CREATE PROBLEMS; MZ - probably need to fix???
!integer, parameter :: taggenericparam=97

integer, parameter :: tagTninf=98
integer, parameter :: tagxnrange=99,tagynrange=104
integer, parameter :: taglx=105,tagxn=106,tagyn=107,tagzn=108,tagdvnx=109


!> VARIABLES REUSED BY ALL WORKERS AND USING MODULES
integer, protected :: myid,lid
!! no external procedure should mess with these (but they need to be able to read them)

integer, private :: ierr=0
!> using procedures need to be able to overwrite this to prevent seg. faults (or something)


!> VARIABLES RELATED TO PROCESS GRID (IF USED)
integer, protected :: lid2,lid3,myid2,myid3


!> Some explanation as the the naming convention used in this module is in order at this point.
!> Generally it is:
!>  <optype>_<send,recv><dims>_<mpi dims>_<optional indicator>

!
!!> THESE INTERFACES OVERLOAD THE MPI GATHER,BROADCAST SUBROUTINES FOR ARRAYS OF DIFFERENT RANKS.
!!> THESE ARE ALSO USEFUL FOR SUBBING IN DIFFERENT SCENARIOS - 1D VS. 2D MPI DIVISIONS ETC.
!interface gather_recv
!  module procedure gather_recv2D_23, gather_recv3D_23, gather_recv4D_23
!end interface gather_recv
!
!interface gather_send
!  module procedure gather_send2D_23, gather_send3D_23, gather_send4D_23
!end interface gather_send
!
!interface bcast_send
!  module procedure bcast_send2D_23, bcast_send3D_23, bcast_send4D_23
!end interface bcast_send
!
!interface bcast_recv
!  module procedure bcast_recv2D_23, bcast_recv3D_23, bcast_recv4D_23
!end interface bcast_recv
!
!interface bcast_send1D_2
!  module procedure bcast_send1D_23_2
!end interface bcast_send1D_2
!interface bcast_recv1D_2
!  module procedure bcast_recv1D_23_2
!end interface bcast_recv1D_2
!
!interface bcast_send1D_3
!  module procedure bcast_send1D_23_3
!end interface bcast_send1D_3
!interface bcast_recv1D_3
!  module procedure bcast_recv1D_23_3
!end interface bcast_recv1D_3
!
!!> THIS ALLOWS EASY SWAPPING OF DIFFERENT ROUTINES FOR 3 VS. 23 DIVISIONS
!interface halo
!  module procedure halo_23
!end interface halo
!interface bcast_send3D_x3i
!  module procedure bcast_send3D_x3i_23
!end interface bcast_send3D_x3i
!interface bcast_recv3D_x3i
!  module procedure bcast_recv3D_x3i_23
!end interface bcast_recv3D_x3i
!
!interface bcast_send3D_x2i
!  module procedure bcast_send3D_x2i_23
!end interface bcast_send3D_x2i
!interface bcast_recv3D_x2i
!  module procedure bcast_recv3D_x2i_23
!end interface bcast_recv3D_x2i
!
!interface bcast_send3D_ghost
!  module procedure bcast_send3D_ghost_23
!end interface bcast_send3D_ghost
!interface bcast_recv3D_ghost
!  module procedure bcast_recv3D_ghost_23
!end interface bcast_recv3D_ghost
!
!interface halo_end
!  module procedure halo_end_23
!end interface halo_end


!> THESE INTERFACES OVERLOAD THE MPI GATHER,BROADCAST SUBROUTINES FOR ARRAYS OF DIFFERENT RANKS.
interface gather_recv
  procedure gather_recv2D_23, gather_recv3D_23, gather_recv4D_23
end interface gather_recv

interface gather_send
  procedure gather_send2D_23, gather_send3D_23, gather_send4D_23
end interface gather_send

interface bcast_send
  procedure bcast_send1D_23, bcast_send2D_23, bcast_send3D_23, bcast_send4D_23
end interface bcast_send

interface bcast_recv
  procedure bcast_recv1D_23, bcast_recv2D_23, bcast_recv3D_23, bcast_recv4D_23
end interface bcast_recv

interface bcast_send1D_2
  module procedure bcast_send1D_23_2
end interface bcast_send1D_2
interface bcast_recv1D_2
  module procedure bcast_recv1D_23_2
end interface bcast_recv1D_2

interface bcast_send1D_3
  module procedure bcast_send1D_23_3
end interface bcast_send1D_3
interface bcast_recv1D_3
  module procedure bcast_recv1D_23_3
end interface bcast_recv1D_3

interface halo
  module procedure halo_23
end interface halo
interface bcast_send3D_x3i
  module procedure bcast_send3D_x3i_23
end interface bcast_send3D_x3i
interface bcast_recv3D_x3i
  module procedure bcast_recv3D_x3i_23
end interface bcast_recv3D_x3i

interface bcast_send3D_x2i
  module procedure bcast_send3D_x2i_23
end interface bcast_send3D_x2i
interface bcast_recv3D_x2i
  module procedure bcast_recv3D_x2i_23
end interface bcast_recv3D_x2i

interface bcast_send3D_ghost
  module procedure bcast_send3D_ghost_23
end interface bcast_send3D_ghost
interface bcast_recv3D_ghost
  module procedure bcast_recv3D_ghost_23
end interface bcast_recv3D_ghost

interface halo_end
  module procedure halo_end_23
end interface halo_end


interface ! mpisend

module subroutine gather_send2D_23(paramtrim,tag)
real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
end subroutine gather_send2D_23

module subroutine gather_send3D_23(paramtrim,tag)
real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
end subroutine gather_send3D_23

module subroutine gather_send4D_23(param,tag)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
end subroutine gather_send4D_23

module subroutine bcast_send1D_23(paramall,tag,param)
real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param
end subroutine bcast_send1D_23

module subroutine bcast_send2D_23(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrim
end subroutine bcast_send2D_23

module subroutine bcast_send3D_23(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim
end subroutine bcast_send3D_23

module subroutine bcast_send4D_23(paramall,tag,param)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
end subroutine bcast_send4D_23

module subroutine bcast_send3D_x3i_23(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim
end subroutine bcast_send3D_x3i_23

module subroutine bcast_send3D_ghost_23(paramall,tag,param)
real(wp), dimension(-1:,-1:,-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
end subroutine bcast_send3D_ghost_23

module subroutine bcast_send3D_x2i_23(paramtrimall,tag,paramtrim)
real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim
end subroutine bcast_send3D_x2i_23

module subroutine bcast_send1D_23_3(paramall,tag,param)
real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param
end subroutine bcast_send1D_23_3

module subroutine bcast_send1D_23_2(paramall,tag,param)
real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param
end subroutine bcast_send1D_23_2

end interface


interface ! mpirecv

module subroutine gather_recv2D_23(paramtrim,tag,paramtrimall)
real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrimall
end subroutine gather_recv2D_23

module subroutine gather_recv3D_23(paramtrim,tag,paramtrimall)
real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrimall
end subroutine gather_recv3D_23

module subroutine gather_recv4D_23(param,tag,paramall)
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: paramall
end subroutine gather_recv4D_23

module subroutine bcast_recv1D_23(param,tag)
real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv1D_23

module subroutine bcast_recv2D_23(paramtrim,tag)
real(wp), dimension(:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv2D_23

module subroutine bcast_recv3D_23(paramtrim,tag)
real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv3D_23

module subroutine bcast_recv4D_23(param,tag)
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv4D_23

module subroutine bcast_recv3D_x3i_23(paramtrim,tag)
real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv3D_x3i_23

module subroutine bcast_recv3D_ghost_23(param,tag)
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv3D_ghost_23

module subroutine bcast_recv1D_old3(param,tag)
real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv1D_old3

module subroutine bcast_recv1D_23_2(param,tag)
real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv1D_23_2

module subroutine bcast_recv1D_23_3(param,tag)
real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag
end subroutine bcast_recv1D_23_3

module subroutine bcast_recv2D_23_3(paramtrim,tag)
real(wp), dimension(:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv2D_23_3

module subroutine bcast_recv3D_x2i_23(paramtrim,tag)
real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag
end subroutine bcast_recv3D_x2i_23

end interface


interface ! mpihalo

module subroutine halo_23(param,lhalo,tag,isperiodic)
  real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
  integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (1 or 2 only)
  integer, intent(in) :: tag
  logical, intent(in) :: isperiodic
end subroutine halo_23

module subroutine halo_end_23(param,paramend,paramtop,tag)
  real(wp), dimension(:,:,:), intent(inout) :: param
  real(wp), dimension(:,:), intent(out) :: paramend
  real(wp), dimension(:,:), intent(out) :: paramtop
  integer, intent(in) :: tag
end subroutine halo_end_23

end interface


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


function grid2ID(i2,i3)

  !------------------------------------------------------------
  !-------COMPUTES A PROCESS ID FROM A LOCATION ON THE PROCESS
  !-------GRID
  !------------------------------------------------------------

  integer, intent(in) :: i2,i3
  integer :: grid2ID

  grid2ID=i3*lid2+i2    !this formulat assumes that the first element is (i2,i3)=(0,0)

end function grid2ID


function ID2grid(ID)

  !------------------------------------------------------------
  !-------COMPUTES GRID LOCATION FROM A PROCESS ID
  !------------------------------------------------------------

  integer, intent(in) :: ID
  integer, dimension(2) :: ID2grid

  ID2grid(2)=ID/lid2                !x3 index into process grid
  ID2grid(1)=ID-ID2grid(2)*lid2     !x2 index into process grid

end function ID2grid


subroutine mpi_manualgrid(lx2all,lx3all,lid2in,lid3in)

integer, intent(in) :: lx2all,lx3all
integer, intent(in) :: lid2in,lid3in

integer, dimension(2) :: inds


if (lx2all/lid2in*lid2in/=lx2all) error stop 'user input grid split in x2 will not work'

if (lx3all/lid3in*lid3in/=lx3all) error stop 'user input grid split in x3 will not work'

if (lid2in*lid3in/=lid) error stop 'total number of processes not commensurate with x2 and x3 split'

lid2=lid2in
lid3=lid3in

!THIS PROCESS' LOCATION ON THE GRID
inds=ID2grid(myid)
myid2=inds(1)
myid3=inds(2)

print *, 'Input process grid is x2 by x3 size (in number of processes):  ',lid2,' by ',lid3
print *, 'Process:  ',myid,' is at location:  ',myid2,myid3,' on the process grid'

end subroutine mpi_manualgrid


subroutine mpigrid(lx2all,lx3all)

!! THIS SUBROUTINE DEFINES A PROCESS GRID, IF REQUIRED

integer, intent(in) :: lx2all,lx3all

integer, dimension(2) :: inds
logical :: x2div,x3div


if (lx3all==1 .or. lx2all==1) then    !this is a 2D simulation so the mpi'd dimenions will get swapped to x3
  lid3=lid         !just divide in x3
  lid2=1
else
  if (lx3all/lid*lid/=lx3all) then
    write (stderr, *) 'lx3all:',lx3all,'  lid:',lid
    error stop 'Grid is not divisible by number of processes (lx3all/lid*lid /= lx3all).'
  end if

  lid2=1
  lid3=lid
  do while( ((lid3/2)*2==lid3) .and. (lid3-lid2>lid3 .or. lid3-lid2>lid2) .and. &
           lx3all/(lid3/2)*(lid3/2)==lx3all .and. lx2all/(lid2*2)*(lid2*2)==lx2all .and. &
           lid3/2>1)
  !! must insure that lx3 is divisible by lid3 and lx2 by lid2 and lid3 must be > 1

    lid3=lid3/2
    lid2=lid2*2
  end do
end if


!FORCE THE CODE TO USE 1D PROCESS GRID
!lid2=1; lid3=lid;


!THIS PROCESS' LOCATION ON THE GRID
inds=ID2grid(myid)
myid2=inds(1)
myid3=inds(2)

print *, 'Proposed process grid is x2 by x3 size (in number of processes):  ',lid2,' by ',lid3
print *, 'Process:  ',myid,' is at location:  ',myid2,myid3,' on the process grid'

end subroutine mpigrid


function slabinds(ID,lx2,lx3)

!! GET THE MIN AND MAX X2,X3 INDICES REFERENCING FULL GRID VARIABLE FOR A GIVEN
!! PROCESS ID

integer, intent(in) :: ID
integer, intent(in) :: lx2,lx3

integer :: i2,i3,i2start,i2fin,i3start,i3fin
integer, dimension(2) :: inds

integer, dimension(4) :: slabinds


  inds=ID2grid(ID)   !find the location on the process grid for this particular process ID
  i2=inds(1)          !need process grid location in order to know where to put the incoming data
  i3=inds(2)
  i3start=i3*lx3+1   !index (3rd dim) in the full grid variable into which the next chunk of data are to be store
  i3fin=i3start+lx3-1
  i2start=i2*lx2+1   !index into 2nd dim of process grid
  i2fin=i2start+lx2-1
  slabinds(1)=i2start
  slabinds(2)=i2fin
  slabinds(3)=i3start
  slabinds(4)=i3fin

end function slabinds


subroutine mpibreakdown()
!! SHUTS DOWN MPI

call mpi_finalize(ierr)

end subroutine mpibreakdown


end module mpimod
