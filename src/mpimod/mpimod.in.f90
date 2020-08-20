module mpimod
!! NOTES:
!! * Need to consider overloading routines as send_ghost and send_noghost so that
!!   it is more clear what the structure of the input arrays should be.
use, intrinsic:: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only : lsp, wp
!! code needs to know how many species are being used.

! use mpi, only: mpi_init, mpi_comm_rank, mpi_comm_size, mpi_comm_world, &
!   mpi_integer,mpi_sum, &
!   mpi_status_size, mpi_status_ignore, MPI_PROC_NULL, &
!   mpi_realprec=>@mpi_realprec@

implicit none !(type, external)

include 'mpif.h'

private
public :: gemini_mpi, myid, myid2, myid3, lid, lid2, lid3, &
  mpi_realprec, mpisetup, mpibreakdown, mpi_manualgrid, mpigrid, id2grid, grid2id, slabinds, &
  bcast_send,  bcast_send1d_2, bcast_send1d_3, bcast_send3d_x2i, bcast_send3d_x3i, bcast_send3d_ghost, &
  bcast_recv, bcast_recv1d_2, bcast_recv1d_3, bcast_recv3d_x2i, bcast_recv3d_x3i, bcast_recv3d_ghost, &
  gather_send, gather_recv, &
  halo, halo_end, &
  mpi_comm_world, mpi_status_ignore, mpi_integer, mpi_sum, &
  test_process_number

!external :: mpi_finalize, mpi_send, mpi_recv, mpi_isend, mpi_irecv, mpi_waitall

integer, parameter :: mpi_realprec = @mpi_realprec@
type :: gemini_mpi_tags

integer :: ns=2, vs1=3, Ts=4
!! root/workers input routines.  also output routines for root/worker
integer :: J1=6,J2=7,J3=8,v2=9,v3=10
!! output root/worker routines, main program, potential_comm
integer :: vs3BC=100,nsBC=101,rhovs1BC=102,rhoesBC=103
!! used in the advection boundary conditions
integer:: vs1BC=1000,vs2BC=1001!,vs3BC=1002
!! used in the compression solution

!> THESE MESSAGES ARE USED IN ELECTRODYNAMICS MODULE
! integer :: E1=11,E2=12,E3=13
integer :: sigP=16,sigH=17,sig0=18,incap=19,v2pol=20,v3pol=21
integer :: DE2Dt=22,DE3Dt=23,flagdirich=24    !unused in present version
! integer :: vn2=25,vn3=26,B1=27                !for passing/gathering full-grid winds

!> IN THE MAIN PROGRAM
integer :: x3=1,dt=5
integer :: x1=27,x2=28

!> IN THE GRID MODULE
integer :: h1=29,h2=30,h3=31
integer :: glat=32,glon=33,alt=34
integer :: lx1=35,lx2=36,lx3=37,lx3all=38
integer :: Bmag=39,inc=40,null=41
integer :: eunit1=42,eunit2=43,eunit3=44,er=45,etheta=46,ephi=47
integer :: r=56,theta=57,phi=58

!> IN THE NEUTRAL MODULE
integer :: lrho=48,lz=49
integer :: dnO=50,dnN2=51,dnO2=52,dTn=53,dvnrho=54,dvnz=55
integer :: ly=69

!> FOR DEALING WITH PRECIPITATION BOUNDARY CONDITIONS MODULE
integer :: llat=59,llon=60,mlat=61,mlon=62,Qp=63,E0p=64

!> FOR DEALING WITH THE ELECTRIC FIELD BOUNDARY CONDITIONS
integer :: E0xp=65,E0yp=66,E0xi=67,E0yi=68

!> FOR DISTRIBUTING PART OF THE ELECTRODYNAMICS CALCULATIONS
integer :: src=69,SigPint2=70,SigPint3=71,SigHint=72,incapint=73,v2electro=74,v3electro=75
integer :: E01=76,E02=77,E03=78,Vminx1=79,Vmaxx1=80

!> THESE ARE USED IN MAGCALC.F90 PROGRAM
integer :: Br=81,Btheta=82, Bphi=83
integer :: dV=84,Jx=85,Jy=86,Rx=87,Ry=88,Rz=89,Rcubed=90,Jz=91

!> FOR COMMUNICATING IF THE GRID DIMENSIONS HAVE BEEN SWAPPED
integer :: swap=92

!> FOR SENDING THE FULL X2 GRID SIZE
integer :: lx2all=93
integer :: x2all=94
integer :: x3all=95

!> AURORAL (S)
integer :: Aur=96

!!> GENERIC PARAMETER (USED BY ADVECTION CODE - HOPEFULLY DOESN'T CREATE PROBLEMS; MZ - probably need to fix???
!integer :: genericparam=97

integer :: Tninf=98
integer :: xnrange=99,ynrange=104
integer :: lx=105,xn=106,yn=107,zn=108,dvnx=109

end type gemini_mpi_tags

type(gemini_mpi_tags), protected :: gemini_mpi

!> A LIST OF TAGS SO THESE DO NOT NEED TO BE EMBEDDED IN EACH SUBROUTINE

!> VARIABLES REUSED BY ALL WORKERS AND USING MODULES
integer, protected :: myid,lid
!! no external procedure should mess with these (but they need to be able to read them)

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


interface ! mpisend.f90
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


interface ! mpirecv.f90
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


interface ! mpihalo.f90
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

integer :: ierr

call mpi_init(ierr)
if (ierr/=0) error stop 'mpimod: mpi_init'
call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
if (ierr/=0) error stop 'mpimod: mpi_comm_rank'
call mpi_comm_size(MPI_COMM_WORLD,lid,ierr)
if (ierr/=0) error stop 'mpimod: mpi_comm_size'

if(myid==0) print *, lid, "MPI processes detected"

!> INITIALIZE, ONLY PARALLELIZING IN X3, GRIDDING FUNCTION MAY CHANGE THIS, IF CALLED.
lid2 = 1
lid3 = lid

end subroutine mpisetup


function slabinds(ID,lx2,lx3)
!! GET THE MIN AND MAX X2,X3 INDICES REFERENCING FULL GRID VARIABLE FOR A GIVEN
!! PROCESS ID

integer, intent(in) :: ID
integer, intent(in) :: lx2,lx3

integer :: i2,i3,i2start,i2fin,i3start,i3fin
integer, dimension(2) :: inds

integer, dimension(4) :: slabinds


inds=ID2grid(ID)
!! find the location on the process grid for this particular process ID
i2=inds(1)
!! need process grid location in order to know where to put the incoming data
i3=inds(2)
i3start=i3*lx3+1
!! index (3rd dim) in the full grid variable into which the next chunk of data are to be store
i3fin=i3start+lx3-1
i2start=i2*lx2+1
!! index into 2nd dim of process grid
i2fin=i2start+lx2-1
slabinds(1)=i2start
slabinds(2)=i2fin
slabinds(3)=i3start
slabinds(4)=i3fin

end function slabinds


subroutine mpi_manualgrid(lx2all,lx3all,lid2in,lid3in)
integer, intent(in) :: lx2all,lx3all, lid2in,lid3in
integer, dimension(2) :: inds

if (lx2all/lid2in*lid2in /= lx2all) error stop 'user input grid split in x2 will not work'
if (lx3all/lid3in*lid3in /= lx3all) error stop 'user input grid split in x3 will not work'
if (lid2in*lid3in /= lid) error stop 'total number of processes not commensurate with x2 and x3 split'

lid2=lid2in
lid3=lid3in

!THIS PROCESS' LOCATION ON THE GRID
inds=ID2grid(myid)
myid2=inds(1)
myid3=inds(2)

end subroutine mpi_manualgrid


subroutine mpigrid(lx2all,lx3all)
!! Automatically determine the PROCESS GRID
!! sets value of lid2,lid3 globally
!! FIXME: should use derived type for these
!! FIXME: improve algorithm to use more CPU cores (be more effective in finding factors for x2 and x3)

integer, intent(in) :: lx2all,lx3all

integer, dimension(2) :: inds

if (lx3all==1) then
  !! 2D simulation, NOT swapped, divide in x3
  lid3 = min(lid, lx2all)
  lid2 = 1
elseif (lx2all==1) then
  !! 2D simulation, SWAP x2 to x3, divide in x3
  lid3 = min(lid, lx3all)
  lid2 = 1
else
  !! 3D simulation
  lid = min(lid, lx3all)
  !! more CPUs than lx3all, reduce used MPI images

  do while(modulo(lx3all, lid) /= 0)
    lid = lid-1
  end do
  !! make number of MPI images a factor of lx3all

  lid2=1
  lid3=lid
  do while( ((lid3/2)*2==lid3) .and. (lid3-lid2>lid3 .or. lid3-lid2>lid2) .and. &
            lx3all/(lid3/2)*(lid3/2)==lx3all .and. lx2all/(lid2*2)*(lid2*2)==lx2all .and. &
            lid3/2>1)
  !! ensure that lx3 is divisible by lid3 and lx2 by lid2 and lid3 must be > 1

    lid3=lid3/2
    lid2=lid2*2
  end do
end if

!> THIS PROCESS' LOCATION ON THE GRID
inds = ID2grid(myid)
myid2 = inds(1)
myid3 = inds(2)

end subroutine mpigrid


integer function grid2id(i2,i3)
!! COMPUTES A PROCESS ID FROM A LOCATION ON THE PROCESS GRID
integer, intent(in) :: i2,i3

grid2ID = i3 * lid2 + i2
!! this formula assumes that the first element is (i2,i3)=(0,0)

end function grid2id


pure function ID2grid(ID)
!! COMPUTES GRID LOCATION FROM A PROCESS ID
integer, dimension(2) :: ID2grid
integer, intent(in) :: ID


ID2grid(2) = ID / lid2
!! x3 index into process grid
ID2grid(1) = ID - ID2grid(2) * lid2
!! x2 index into process grid

end function ID2grid


integer function mpibreakdown() result(ierr)
!! SHUTS DOWN MPI

call mpi_finalize(ierr)

end function mpibreakdown


subroutine test_process_number(N, lx2all, lx3all, rx2, rx3)
!! this is only for testing. Due to lid being protected, has to be in this module
!! FIXME: make type for MPI parameters.

integer, intent(in) :: N(:), rx2(:), rx3(:), lx2all, lx3all
integer :: i

do i = 1,size(N)
  lid = N(i)
  call mpigrid(lx2all,lx3all)
  if (lid2 /= rx2(i) .or. lid3 /= rx3(i)) then
    write(stderr,'(A,5I4)') 'failed: lx2all,lx3all,lid,N:',lx2all,lx3all,lid,N(i)
    write(stderr,*) 'expected lid2,lid3', rx2(i), rx3(i), 'but got:',lid2,lid3
    error stop
  end if
end do

end subroutine test_process_number


end module mpimod
