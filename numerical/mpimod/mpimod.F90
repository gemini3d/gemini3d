module mpimod
!! NOTES:
!! * Need to consider overloading routines as send_ghost and send_noghost so that
!!   it is more clear what the structure of the input arrays should be. 

use phys_consts, only : lsp, wp 
!! code needs to know how many species are being used.

use mpi, only: mpi_init, mpi_comm_world, &
               mpi_integer,mpi_sum, &
               mpi_status_size, mpi_status_ignore, MPI_PROC_NULL, &
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

!> FOR SENDING THE FULL X2 GRID SIZE
integer, parameter :: taglx2all=93
integer, parameter :: tagx2all=94
integer, parameter :: tagx3all=95

!> AURORAL TAG(S)
integer, parameter :: tagAur=96

!!> GENERIC PARAMETER (USED BY ADVECTION CODE - HOPEFULLY DOESN'T CREATE PROBLEMS; MZ - probably need to fix???
!integer, parameter :: taggenericparam=97

integer, parameter :: tagTninf=98

!> VARIABLES REUSED BY ALL WORKERS AND USING MODULES
integer, protected :: myid,lid    
!! no external procedure should mess with these (but they need to be able to read them)

integer :: ierr
!> using procedures need to be able to overwrite this to prevent seg. faults (or something)


!> VARIABLES RELATED TO PROCESS GRID (IF USED)
integer, protected :: lid2,lid3,myid2,myid3


!> Some explanation as the the naming convention used in this module is in order at this point.
!> Generally it is:
!>  <optype>_<send,recv><dims>_<mpi dims>_<optional indicator>

!> THESE INTERFACES OVERLOAD THE MPI GATHER,BROADCAST SUBROUTINES FOR ARRAYS OF DIFFERENT RANKS.
!> THESE ARE ALSO USEFUL FOR SUBBING IN DIFFERENT SCENARIOS - 1D VS. 2D MPI DIVISIONS ETC. 
interface gather_recv
  module procedure gather_recv2D_23, gather_recv3D_23, gather_recv4D_23
end interface gather_recv

interface gather_send
  module procedure gather_send2D_23, gather_send3D_23, gather_send4D_23
end interface gather_send

interface bcast_send
  module procedure bcast_send2D_23, bcast_send3D_23, bcast_send4D_23
end interface bcast_send

interface bcast_recv
  module procedure bcast_recv2D_23, bcast_recv3D_23, bcast_recv4D_23
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

!> THIS ALLOWS EASY SWAPPING OF DIFFERENT ROUTINES FOR 3 VS. 23 DIVISIONS
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
  module procedure halo_end_3
end interface halo_end

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


subroutine mpigrid(lx2all,lx3all)

!! THIS SUBROUTINE DEFINES A PROCESS GRID, IF REQUIRED 

integer, intent(in) :: lx2all,lx3all

integer, dimension(2) :: inds

if (lx3all==1 .or. lx2all==1) then    !this is a 2D simulation so the mpi'd dimenions will get swapped to x3
  lid3=lid         !just divide in x3
  lid2=1
else
  if (lx3all/lid*lid/=lx3all) then
    error stop '!!!Grid is not divisible by number of processes - please generate a new one  &
                 & and try again or try a different number of processes...'
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


subroutine halo_3(param,lhalo,tag)
!! GENERIC HALOING ROUTINE FOR FILLING GHOST CELLS.  CAN
!! BE USED TO SET BOUNDARY CONDITIONS OR PREPARE ARRAYS
!! FOR FINITE DIFFERENCING, ETC.  OBVIOUSLY ARRAYS INCLUDE
!! GHOST CELLS.  ARRAYS SHOULD HAVE SPECIES DIMENSION
!! REMOVED BEFORE PASSAGE INTO THIS SUBROUTINE.  THIS
!! ROUTINE WORKS FOR A 3D ARRAY, AND BY DEFAULT ENFORCES
!! PERIODIC BOUNDARY CONDITIONS.  IF APERIODIC CONDITIONS
!! ARE NEEDED YOU MUST OVERWRITE THE GLOBAL BOUNDARIES AFTER
!! YOU CALL HALO.
!! 
!! THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS

real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (1 or 2 only)
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,ihalo
integer :: idleft,idright

integer, dimension(4) :: requests
integer, dimension(MPI_STATUS_SIZE,4) :: statuses
integer :: tmpreq

real(wp) :: tstart,tfin

lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> IDENTIFY MY NEIGHBORS
idleft=myid-1
idright=myid+1


if (.not. (myid==0 .and. idright==lid)) then    !if root is the only worker, do nothing...
  
  !> SCREEN FOR GLOBAL BOUNDARIES, ASSUME PERIODIC.
  !> MUST BE OVERWRITTEN LATER IF YOU ARE USING ANOTHER TYPE OF BOUNDARY
  if (idleft==-1) then
    idleft=lid-1
  !      idleft=MPI_PROC_NULL    !if you wanted to default to aperiodic you could do this...
  end if
  if (idright==lid) then
    idright=0
  !      idright=MPI_PROC_NULL
  end if
  
  call mpi_isend(param(:,:,1:lhalo),(lx1+4)*(lx2+4)*lhalo,mpi_realprec,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(1)=tmpreq
  call mpi_isend(param(:,:,lx3+1-lhalo:lx3),(lx1+4)*(lx2+4)*lhalo,mpi_realprec, &
                    idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(2)=tmpreq
  call mpi_irecv(param(:,:,lx3+1:lx3+lhalo),(lx1+4)*(lx2+4)*lhalo,mpi_realprec,idright, &
                    tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(3)=tmpreq
  call mpi_irecv(param(:,:,1-lhalo:0),(lx1+4)*(lx2+4)*lhalo,mpi_realprec,idleft, &
                          tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(4)=tmpreq
  
  call mpi_waitall(4,requests,statuses,ierr)
  
end if

end subroutine halo_3


subroutine halo_end_3(param,paramend,tag)

!MZ - needs to be updated to accommodate x2/x3 division...

!! GENERIC HALOING ROUTINE WHICH PASSES THE BEGINNING OF THE
!! SLAB TO ITS LEFTWARD (IN X3) NEIGHBOR SO THAT X3 INTEGRATIONS
!! CAN BE DONE PROPERLY.  PRESENTLY THIS IS JUST USED IN MAGCALC
!! 
!! THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS.  IT
!! ALSO ASSUMES THAT THE ARRAYS INVOLVED DO HAVE GHOST CELLS

real(wp), dimension(:,:,:), intent(inout) :: param
real(wp), dimension(:,:), intent(out) :: paramend
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,ihalo
integer :: idleft,idright

integer, dimension(2) :: requests
integer, dimension(MPI_STATUS_SIZE,4) :: statuses
integer :: tmpreq

real(wp) :: tstart,tfin

lx1=size(param,1)
lx2=size(param,2)
lx3=size(param,3)


!IDENTIFY MY NEIGHBORS, I NEED TO GET DATA FROM BEGINNING OF RIGHT (FOR MY
!END) AND SEND MY BEGINNING DATA TO THE LEFT (FOR THEIR END)
idleft=myid-1
idright=myid+1


if (.not. (myid==0 .and. idright==lid)) then    !if root is the only worker, do nothing...
  
  !SCREEN FOR GLOBAL BOUNDARIES, ASSUME PERIODIC (MUST BE OVERWRITTEN LATER IF
  !YOU ARE USING ANOTHER TYPE OF BOUNDARY
  if (idleft==-1) then
    idleft=lid-1
  end if
  if (idright==lid) then
    idright=0
  end if
  
  
  call mpi_isend(param(:,:,1),lx1*lx2,mpi_realprec,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(1)=tmpreq
  call mpi_irecv(paramend,lx1*lx2,mpi_realprec,idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
  requests(2)=tmpreq
  
  call mpi_waitall(2,requests,statuses,ierr)
  
  if (myid==lid-1) paramend=0d0    !zero out the data at the end of the grid
  
end if

end subroutine halo_end_3


!   LEAVE THIS IN FOR FUTURE DEV. - HALOS AN ARRAY THAT IS SPLIT ALONG THE 2 AND
!   3 RANKS
subroutine halo_23(param,lhalo,tag,isperiodic)

  !------------------------------------------------------------
  !-------GENERIC HALOING ROUTINE FOR FILLING GHOST CELLS.  CAN
  !-------BE USED TO SET BOUNDARY CONDITIONS OR PREPARE ARRAYS
  !-------FOR FINITE DIFFERENCING, ETC.  OBVIOUSLY ARRAYS INCLUDE
  !-------GHOST CELLS.  ARRAYS SHOULD HAVE SPECIES DIMENSION
  !-------REMOVED BEFORE PASSAGE INTO THIS SUBROUTINE.  THIS
  !-------ROUTINE WORKS FOR A 3D ARRAY.
  !-------
  !-------THIS SUBROUTINE SHOULD BE CONSIDERED DISTINCT IN FUNCTIONALITY
  !-------FROM THE SUBROUTINE WHICH COMPUTES BOUNDARY CONDITIONS FOR
  !-------THE GLOBAL GRID.  I.E., IF CALLED ON A WORKER THAT ABUTS THE
  !-------GLOBAL BOUNDARY IT WILL DO *NOTHING*.
  !-------
  !-------THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS
  !-------
  !-------THIS VERSION ALSO ASSUMES A PROCESS GRID HAS BEEN DEFINED
  !-------AND THAT PASSING NEEDS TO BE DONE IN X2 AND X3
  !------------------------------------------------------------

  real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
  integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (probably 1 or 2)
  integer, intent(in) :: tag
  logical, intent(in) :: isperiodic

  integer :: lx1,lx2,lx3,ihalo
  integer :: idleft,idright,idup,iddown
  integer :: i2,i3

  integer, dimension(4) :: requests
  integer, dimension(MPI_STATUS_SIZE,4) :: statuses
  integer :: tmpreq

  real(wp), allocatable, dimension(:,:,:) :: buffer31,buffer32,buffer33,buffer34
  real(wp), allocatable, dimension(:,:,:) :: buffer21,buffer22,buffer23,buffer24
  real(wp) :: tstart,tfin

  logical :: x2begin,x3begin,x2end,x3end     !to store info about whether we are the first or last process in a direction


  !COMPUTE SIZES, JUST IN CASE.  
  lx1=size(param,1)-4
  lx2=size(param,2)-4
  lx3=size(param,3)-4


  !IDENTIFY MY NEIGHBORS IN X3
  x3begin=.false.
  x3end=.false.

  i3=myid3-1
  i2=myid2
  if (i3==-1) then          !global boundary to my left, assume periodic
    i3=lid3-1               !lid3-1 is the last process in x3 on the process grid
    x3begin=.true.
  end if
  idleft=grid2ID(i2,i3)
  if (x3begin .and. .not.(isperiodic)) then     !we are flagged as not wanting periodic boundaries so do nothing (overwrite idleft to send to NULL process
    idleft=MPI_PROC_NULL
  end if


  i3=myid3+1
  i2=myid2
  if (i3==lid3) then        !global boundary to my right, assume periodic
    i3=0
    x3end=.true.
  end if
  idright=grid2ID(i2,i3)    !convert the location on process grid into a flat processed ID, The process grid is 
                            !visualized as lid2,lid3 in terms of index order (e.g. the i2 index cycles more quickly
  if (x3end .and. .not.(isperiodic)) then
    idright=MPI_PROC_NULL
  end if

  !IDENTIFY MY NEIGHBORING PROCESSES IN X2
  x2begin=.false.
  x2end=.false.

  i3=myid3
  i2=myid2-1
  if (i2==-1) then       !global boundary downward, assume periodic
    i2=lid2-1
    x2begin=.true.
  end if
  iddown=grid2ID(i2,i3)
  if (x2begin) then      !never assume periodic in the x2-direction
    iddown=MPI_PROC_NULL
  end if

  i3=myid3
  i2=myid2+1
  if (i2==lid2) then     !global boundary upward, assume periodic
    i2=0
    x2end=.true.
  end if
  idup=grid2ID(i2,i3)    !convert to process ID
  if (x2end) then
    idup=MPI_PROC_NULL
  end if

!  !some debug output
!  print*, 'Computing neighbors for ID:  ',myid,' at location on the process grid:  ',myid2,myid3
!  print*, iddown,idup,idleft,idright


  !ALLOCATE SPACE TO BUFFER MESSAGES (SINCE USING ASYNCHRONOUS MPI COMMANDS), THESE HAVE TO BE ALLOCATED SINCE WE
  !DON'T KNOW A PRIORI HOW MANY CELLS TO HALO.  ALSO NOTE THAT ONLY THE X1-DIRECTION HAS GHOST CELLS FOR THESE.
  allocate(buffer31(-1:lx1+2,1:lx2,lhalo),buffer32(-1:lx1+2,1:lx2,lhalo),buffer33(-1:lx1+2,1:lx2,lhalo), &
           buffer34(-1:lx1+2,1:lx2,lhalo))
  allocate(buffer21(-1:lx1+2,lhalo,1:lx3),buffer22(-1:lx1+2,lhalo,1:lx3),buffer23(-1:lx1+2,lhalo,1:lx3), & 
           buffer24(-1:lx1+2,lhalo,1:lx3))


  !EXCHANGE MESSAGES IN THE X3-DIRECTION OF THE PROCESS GRID
  if (.not. (x3begin .and. x3end)) then        !make sure we actually need to pass in this direction, viz. we aren't both the beginning and thend
    buffer31=param(-1:lx1+2,1:lx2,1:lhalo)     !x1 ghost cells to be included
    call mpi_isend(buffer31,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(1)=tmpreq
  
    buffer32=param(-1:lx1+2,1:lx2,lx3+1-lhalo:lx3)
    call mpi_isend(buffer32,(lx1+4)*(lx2)*lhalo,mpi_realprec, &
                      idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(2)=tmpreq
  
    call mpi_irecv(buffer33,(lx1+4)*(lx2)*lhalo,mpi_realprec,idright, &
                      tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(3)=tmpreq
  
    call mpi_irecv(buffer34,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft, &
                            tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(4)=tmpreq
  
    call mpi_waitall(4,requests,statuses,ierr)

    if (idright/=MPI_PROC_NULL)  then    !only overwrite the cells if we didn't do a null receive
      param(-1:lx1+2,1:lx2,lx3+1:lx3+lhalo)=buffer33    !can't copy out buffers until we know the messages have been received
    end if
    if (idleft/=MPI_PROC_NULL) then
      param(-1:lx1+2,1:lx2,1-lhalo:0)=buffer34
    end if
  end if

  !EXCHANGE MESSAGES IN THE X2 DIRECTION OF THE PROCESS GRID
  if (.not. (x2begin .and. x2end)) then
    buffer21=param(-1:lx1+2,1:lhalo,1:lx3)
    call mpi_isend(buffer21,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(1)=tmpreq
  
    buffer22=param(-1:lx1+2,lx2+1-lhalo:lx2,1:lx3)
    call mpi_isend(buffer22,(lx1+4)*(lx3)*lhalo,mpi_realprec, &
                      idup,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(2)=tmpreq
  
    call mpi_irecv(buffer23,(lx1+4)*(lx3)*lhalo,mpi_realprec,idup,&
                      tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(3)=tmpreq
  
    call mpi_irecv(buffer24,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown, &
                            tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(4)=tmpreq
    
    call mpi_waitall(4,requests,statuses,ierr)

    if (idup/=MPI_PROC_NULL) then
      param(-1:lx1+2,lx2+1:lx2+lhalo,1:lx3)=buffer23    !clear to copy out buffers
    end if
    if (iddown/=MPI_PROC_NULL) then
      param(-1:lx1+2,1-lhalo:0,1:lx3)=buffer24
    end if
  end if

  !CLEAR OUT BUFFER VARIABLES
  deallocate(buffer31,buffer32,buffer33,buffer34)
  deallocate(buffer21,buffer22,buffer23,buffer24) 

end subroutine halo_23


subroutine halo_end_23(param,paramend,paramtop,tag)

  !! GENERIC HALOING ROUTINE WHICH PASSES THE BEGINNING OF THE
  !! SLAB TO ITS LEFTWARD (IN X3) NEIGHBOR SO THAT X3 INTEGRATIONS
  !! CAN BE DONE PROPERLY.  PRESENTLY THIS IS JUST USED IN MAGCALC
  !! 
  !! THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS.  IT
  !! ALSO ASSUMES THAT THE ARRAYS INVOLVED DO HAVE GHOST CELLS
  
  real(wp), dimension(:,:,:), intent(inout) :: param
  real(wp), dimension(:,:), intent(out) :: paramend
  real(wp), dimension(:,:), intent(out) :: paramtop
  integer, intent(in) :: tag
  
  integer :: lx1,lx2,lx3,ihalo
  integer :: idleft,idright,iddown,idup
  integer :: i2,i3
  
  integer, dimension(2) :: requests
  integer, dimension(MPI_STATUS_SIZE,4) :: statuses
  integer :: tmpreq
  
  real(wp) :: tstart,tfin
  logical :: x2begin,x2end,x3begin,x3end
  real(wp), dimension(:,:), allocatable :: buffer
 
 
  lx1=size(param,1)
  lx2=size(param,2)
  lx3=size(param,3)
  

  !IDENTIFY MY NEIGHBORS, I NEED TO GET DATA FROM BEGINNING OF RIGHT (FOR MY
  !END) AND SEND MY BEGINNING DATA TO THE LEFT (FOR THEIR END)
  !IDENTIFY MY NEIGHBORS IN X3
  x3begin=.false.
  x3end=.false.

  i3=myid3-1
  i2=myid2
  if (i3==-1) then          !global boundary to my left, assume periodic
    i3=lid3-1               !lid3-1 is the last process in x3 on the process grid
    x3begin=.true.
  end if
  idleft=grid2ID(i2,i3)
  if (x3begin) then     !we are flagged as not wanting periodic boundaries so do nothing (overwrite idleft to send to NULL process
    idleft=MPI_PROC_NULL
  end if


  i3=myid3+1
  i2=myid2
  if (i3==lid3) then        !global boundary to my right, assume periodic
    i3=0
    x3end=.true.
  end if
  idright=grid2ID(i2,i3)    !convert the location on process grid into a flat processed ID, The process grid is 
                            !visualized as lid2,lid3 in terms of index order (e.g. the i2 index cycles more quickly
  if (x3end) then
    idright=MPI_PROC_NULL
  end if

  !IDENTIFY MY NEIGHBORING PROCESSES IN X2
  x2begin=.false.
  x2end=.false.

  i3=myid3
  i2=myid2-1
  if (i2==-1) then       !global boundary downward, assume periodic
    i2=lid2-1
    x2begin=.true.
  end if
  iddown=grid2ID(i2,i3)
  if (x2begin) then      !never assume periodic in the x2-direction
    iddown=MPI_PROC_NULL
  end if

  i3=myid3
  i2=myid2+1
  if (i2==lid2) then     !global boundary upward, assume periodic
    i2=0
    x2end=.true.
  end if
  idup=grid2ID(i2,i3)    !convert to process ID
  if (x2end) then
    idup=MPI_PROC_NULL
  end if


  !PASS DATA IN X3 DIRECTION
  if (.not. (x3begin .and. x3end)) then        !make sure we actually need to pass in this direction, viz. we aren't both the beginning and thend  
    call mpi_isend(param(:,:,1),lx1*lx2,mpi_realprec,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(1)=tmpreq
    call mpi_irecv(paramend,lx1*lx2,mpi_realprec,idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(2)=tmpreq
  
    call mpi_waitall(2,requests,statuses,ierr)
  end if
  

  !PASS DATA IN X3 DIRECTION
  if (.not. (x2begin .and. x2end)) then
    allocate(buffer(lx1,lx3))
    buffer=param(:,1,:)
    call mpi_isend(buffer,lx1*lx3,mpi_realprec,iddown,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(1)=tmpreq
    call mpi_irecv(paramtop,lx1*lx3,mpi_realprec,idup,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(2)=tmpreq
    deallocate(buffer)

    call mpi_waitall(2,requests,statuses,ierr)
  end if


  !ZERO OUT THE ENDS (DO NOT ADD DATA PAST EDGE OF THE GRID
  if (myid2==lid2-1) paramtop=0d0    !add nothing on the end...
  if (myid3==lid3-1) paramend=0d0    !zero out the data at the end of the grid
  
end subroutine halo_end_23


subroutine gather_recv2D_3(paramtrim,tag,paramtrimall)
!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!! 
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!! 
!! THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE ANY GHOST CELLS!!!!

real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrimall

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

end subroutine gather_recv2D_3


subroutine gather_recv2D_23(paramtrim,tag,paramtrimall)
!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!! 
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!! 
!! THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE ANY GHOST CELLS!!!!

real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrimall

integer :: lx1,lx2,lx3,lsp,lx2all,lx3all
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2)) :: paramtmp

lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx3=size(paramtrim,2)


!PATCH DATA TOGETHER FOR OUTPUT STARTING WITH ROOT'S SLAB
paramtrimall(1:lx2,1:lx3)=paramtrim   !copy root's data into full-grid array

do iid=1,lid-1
  call mpi_recv(paramtmp,lx2*lx3, &
                mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  inds=slabinds(iid,lx2,lx3)
  paramtrimall(inds(1):inds(2),inds(3):inds(4))=paramtmp    !note the exclusion of the ghost cells
end do

end subroutine gather_recv2D_23


subroutine gather_recv3D_3(paramtrim,tag,paramtrimall)
!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!! 
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!! 
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! ANY GHOST CELLS!!!!

real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrimall

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

end subroutine gather_recv3D_3


subroutine gather_recv3D_23(paramtrim,tag,paramtrimall)

!! THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!! A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!! OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!! 
!! THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
!! 
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! ANY GHOST CELLS!!!!
!! THIS VERION ALSO WORKS ON A PROCESS GRID

real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrimall

integer :: lx1,lx2,lx3,lx2all,lx3all
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp   !buffer space for mpi receive, includes only x1 ghost cells


lx1=size(paramtrim,1)
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!Originally the outer loop was over worker number, which cycles the 3rd dimension
!slower than 4th (backward from what would be most efficient memory access pattern)
!Since the gathering operation is root-limited probably, I'm guessing it's better
!to give root an efficient memory access pattern here, but I haven't tested this
!theory.
paramtrimall(:,1:lx2,1:lx3)=paramtrim(:,1:lx2,1:lx3)    !store root's piece of data
do iid=1,lid-1        !must loop over all processes in the grid, don't enter loop if only root is present
  call mpi_recv(paramtmp,lx1*lx2*lx3, &          !note no ghost cells!!!
                mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)    !recieve chunk of data into buffer
  inds=slabinds(iid,lx2,lx3)
  paramtrimall(1:lx1,inds(1):inds(2),inds(3):inds(4))=paramtmp    !note the exclusion of the ghost cells
end do

end subroutine gather_recv3D_23


subroutine gather_recv4D_3(param,tag,paramall)

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

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: paramall

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

end subroutine gather_recv4D_3


subroutine gather_recv4D_23(param,tag,paramall)

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

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: paramall

integer :: lx1,lx2,lx3,isp,lx2all,lx3all
integer :: iid
integer, dimension(4) :: inds
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp   !buffer space for mpi receive, includes only x1 ghost cells


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4

!Originally the outer loop was over worker number, which cycles the 3rd dimension
!slower than 4th (backward from what would be most efficient memory access pattern)
!Since the gathering operation is root-limited probably, I'm guessing it's better
!to give root an efficient memory access pattern here, but I haven't tested this
!theory.
do isp=1,lsp
  paramall(-1:lx1+2,1:lx2,1:lx3,isp)=param(-1:lx1+2,1:lx2,1:lx3,isp)    !root records his own piece of the grid into full grid variable

  do iid=1,lid-1        !must loop over all processes in the grid, don't enter loop if only root is present
    call mpi_recv(paramtmp,(lx1+4)*lx2*lx3, &
                  mpi_realprec,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)    !recieve chunk of data into buffer
    inds=slabinds(iid,lx2,lx3)
    paramall(-1:lx1+2,inds(1):inds(2),inds(3):inds(4),isp)=paramtmp(-1:lx1+2,1:lx2,1:lx3)    !note the inclusion of x1 ghost cells
  end do
end do

end subroutine gather_recv4D_23


subroutine gather_send2D_3(paramtrim,tag)

!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
!-------ANY GHOST CELLS!!!!
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag

integer :: lx2,lx3


lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx3=size(paramtrim,2)

call mpi_send(paramtrim,lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)

end subroutine gather_send2D_3


subroutine gather_send2D_23(paramtrim,tag)

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

real(wp), dimension(:,:), intent(in) :: paramtrim
integer, intent(in) :: tag

integer :: lx2,lx3


lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx3=size(paramtrim,2)

call mpi_send(paramtrim,lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)

end subroutine gather_send2D_23


subroutine gather_send3D_3(paramtrim,tag)

!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------ANY GHOST CELLS!!!!
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)

call mpi_send(paramtrim,lx1*lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)

end subroutine gather_send3D_3


subroutine gather_send3D_23(paramtrim,tag)

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

real(wp), dimension(:,:,:), intent(in) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)

call mpi_send(paramtrim,lx1*lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)

end subroutine gather_send3D_23


subroutine gather_send4D_3(param,tag)

!------------------------------------------------------------
!-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
!-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
!-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO GATHER
!-------
!-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,isp


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4

do isp=1,lsp
  call mpi_send(param(:,:,1:lx3,isp),(lx1+4)*(lx2+4)*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)
end do

end subroutine gather_send4D_3


subroutine gather_send4D_23(param,tag)

!------------------------------------------------------------
!-------SENDS 4D DATA ON A 2D PROCESS GRID TO ROOT.
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,isp
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4

do isp=1,lsp
  paramtmp=param(-1:lx1+2,1:lx2,1:lx3,isp)
  call mpi_send(paramtmp,(lx1+4)*lx2*lx3,mpi_realprec,0,tag,MPI_COMM_WORLD,ierr)
end do

end subroutine gather_send4D_23


subroutine bcast_send1D_3(paramall,tag,param)

!------------------------------------------------------------
!-------BROADCASTS MPI DIMENSION VARIABLES TO WORKERS.  NOTE THAT
!-------WE'VE ELECTED TO NOT USE THE GENERAL BROADCAST ROUTINES FOR
!-------SINCE THESE OPERATIONS REQUIRE A LOT OF SPECIAL CASING FOR
!-------THE SIZES OF THE VARIABLES TO BE SENT
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 1D ARRAYS
!------------------------------------------------------------

real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param

integer :: lx,lxall     !local sizes
integer :: iid,islstart,islfin


lxall=size(paramall,1)-4
lx=size(param,1)-4


do iid=1,lid-1
  islstart=iid*lx+1
  islfin=islstart+lx-1

  call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do
param=paramall(-1:lx+2)

end subroutine bcast_send1D_3


subroutine bcast_send1D_23_2(paramall,tag,param)

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

real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param

integer :: lx,lxall     !local sizes
integer :: iid,islstart,islfin
integer, dimension(2) :: indsgrid

lxall=size(paramall,1)-4
lx=size(param,1)-4


do iid=1,lid-1
  indsgrid=ID2grid(iid)
  
!  islstart=iid*lx+1
!  islfin=islstart+lx-1
  islstart=indsgrid(1)*lx+1
  islfin=islstart+lx-1

  call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do
param=paramall(-1:lx+2)

end subroutine bcast_send1D_23_2


subroutine bcast_send1D_23_3(paramall,tag,param)

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

real(wp), dimension(-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:), intent(out) :: param

integer :: lx,lxall     !local sizes
integer :: iid,islstart,islfin
integer, dimension(2) :: indsgrid

lxall=size(paramall,1)-4
lx=size(param,1)-4


do iid=1,lid-1
  indsgrid=ID2grid(iid)    !compute my location on the process grid

!  islstart=iid*lx+1
!  islfin=islstart+lx-1
  islstart=indsgrid(2)*lx+1    !send piece of grid that corresponds to my x3 position
  islfin=islstart+lx-1

  call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do
param=paramall(-1:lx+2)

end subroutine bcast_send1D_23_3


subroutine bcast_send2D_3(paramtrimall,tag,paramtrim)

!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrim

integer :: lx2,lx3
integer :: iid,islstart,islfin


lx2=size(paramtrim,1)    !assume this is an array which has been 'flattened' along the 1-dimension
lx3=size(paramtrim,2)


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_send(paramtrimall(:,islstart:islfin),lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,1:lx3)

end subroutine bcast_send2D_3


subroutine bcast_send2D_23(paramtrimall,tag,paramtrim)

!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:), intent(out) :: paramtrim

integer :: lx2,lx3
integer :: iid,islstart,islfin
integer, dimension(4) :: inds

real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2)) :: paramtmp


lx2=size(paramtrim,1)    !assume this is an array which has been 'flattened' along the 1-dimension
lx3=size(paramtrim,2)


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
!  islstart=iid*lx3+1
!  islfin=islstart+lx3-1
  inds=slabinds(iid,lx2,lx3)

  paramtmp=paramtrimall(inds(1):inds(2),inds(3):inds(4))
  call mpi_send(paramtmp,lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim(1:lx2,1:lx3)=paramtrimall(1:lx2,1:lx3)

end subroutine bcast_send2D_23


subroutine bcast_send3D_3(paramtrimall,tag,paramtrim)

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

real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_send(paramtrimall(:,:,islstart:islfin),lx1*lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!> ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,:,1:lx3)

end subroutine bcast_send3D_3


subroutine bcast_send3D_23(paramtrimall,tag,paramtrim)

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

real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(1:lx1,inds(1):inds(2),inds(3):inds(4))
  call mpi_send(paramtmp,lx1*lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!> ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(1:lx1,1:lx2,1:lx3)

end subroutine bcast_send3D_23


subroutine bcast_send3D_x3i_3(paramtrimall,tag,paramtrim)

!------------------------------------------------------------
!-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!-------
!-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!-------GHOST CELLS, BUT ARE X3 INTERFACE QUANITTIES
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin


lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)-1    !note that we are interpreting input as an x3i quantity meaning that it has size lx3+1


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_send(paramtrimall(:,:,islstart:islfin+1),lx1*lx2*(lx3+1), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)     !note the +1 since thes are interfact quantities (and need to overlap b/t workers)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,:,1:lx3+1)

end subroutine bcast_send3D_x3i_3


subroutine bcast_send3D_x3i_23(paramtrimall,tag,paramtrim)

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

real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp    !has size lx3+1 due to input having that size

lx1=size(paramtrim,1)      !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)-1    !note that we are interpreting input as an x3i quantity meaning that it has size lx3+1


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(:,inds(1):inds(2),inds(3):inds(4)+1)     !+1 since this is an x3 interface quantity
  call mpi_send(paramtmp,lx1*lx2*(lx3+1), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)         !note the +1 since thes are interface quantities (and need to overlap b/t workers)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,1:lx2,1:lx3+1)

end subroutine bcast_send3D_x3i_23


subroutine bcast_send3D_x2i_23(paramtrimall,tag,paramtrim)

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

real(wp), dimension(:,:,:), intent(in) :: paramtrimall
integer, intent(in) :: tag
real(wp), dimension(:,:,:), intent(out) :: paramtrim

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin
integer, dimension(4) :: inds
real(wp), dimension(1:size(paramtrim,1),1:size(paramtrim,2),1:size(paramtrim,3)) :: paramtmp    !has size lx3+1 due to input having that size

lx1=size(paramtrim,1)      !note here that paramtrim does not have ghost cells
lx2=size(paramtrim,2)-1    !for this routine the interface is in the x2-direciton
lx3=size(paramtrim,3)


!ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramtrimall(:,inds(1):inds(2)+1,inds(3):inds(4))     !+1 since this is an x3 interface quantity
  call mpi_send(paramtmp,lx1*(lx2+1)*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)         !note the +1 since thes are interface quantities (and need to overlap b/t workers)
end do


!ROOT TAKES A SLAB OF DATA
paramtrim=paramtrimall(:,1:lx2+1,1:lx3)

end subroutine bcast_send3D_x2i_23


subroutine bcast_send3D_ghost_3(paramall,tag,param)
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!!ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH INCLUDE GHOST CELLS

real(wp), dimension(-1:,-1:,-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin

!> note here that param has ghost cells
lx1=size(param,1)-4    
lx2=size(param,2)-4
lx3=size(param,3)-4


!> ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  islstart=iid*lx3+1
  islfin=islstart+lx3-1

  call mpi_send(paramall(:,:,islstart-2:islfin+2),(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!> ROOT TAKES A SLAB OF DATA
param=paramall(:,:,-1:lx3+2)

end subroutine bcast_send3D_ghost_3


subroutine bcast_send3D_ghost_23(paramall,tag,param)
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!!ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH INCLUDE GHOST CELLS

real(wp), dimension(-1:,-1:,-1:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:), intent(out) :: param

integer :: lx1,lx2,lx3
integer :: iid,islstart,islfin
integer, dimension(4) :: inds
real(wp), dimension(-1:size(param,1)-2,-1:size(param,2)-2,-1:size(param,3)-2) :: paramtmp

!> note here that param has ghost cells
lx1=size(param,1)-4   
lx2=size(param,2)-4
lx3=size(param,3)-4


!> ROOT BROADCASTS IC DATA TO WORKERS
do iid=1,lid-1
  inds=slabinds(iid,lx2,lx3)
  paramtmp=paramall(:,inds(1)-2:inds(2)+2,inds(3)-2:inds(4)+2)
  call mpi_send(paramtmp,(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
end do


!> ROOT TAKES A SLAB OF DATA
param=paramall(:,-1:lx2+2,-1:lx3+2)

end subroutine bcast_send3D_ghost_23


subroutine bcast_send4D_3(paramall,tag,param)
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!! ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param

integer :: lx1,lx2,lx3,isp
integer :: iid,islstart,islfin


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> ROOT BROADCASTS IC DATA TO WORKERS
do isp=1,lsp
  param(:,:,:,isp)=paramall(:,:,-1:lx3+2,isp)
    !! roots part of the data

  do iid=1,lid-1
    islstart=iid*lx3+1
    islfin=islstart+lx3-1

    call mpi_send(paramall(:,:,islstart-2:islfin+2,isp),(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
  end do
end do

end subroutine bcast_send4D_3


subroutine bcast_send4D_23(paramall,tag,param)
!! THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
!! ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
!!
!! SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
integer, intent(in) :: tag
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param

integer :: lx1,lx2,lx3,isp
integer :: iid,islstart,islfin
integer, dimension(4) :: inds
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp


lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> ROOT BROADCASTS IC DATA TO WORKERS
do isp=1,lsp
  param(:,1:lx2,1:lx3,isp)=paramall(:,1:lx2,1:lx3,isp)    ! roots part of the data

  do iid=1,lid-1
    inds=slabinds(iid,lx2,lx3)
    paramtmp=paramall(-1:lx1+2,inds(1):inds(2),inds(3):inds(4),isp)
    call mpi_send(paramtmp,(lx1+4)*lx2*lx3, &
               mpi_realprec,iid,tag,MPI_COMM_WORLD,ierr)
  end do
end do

end subroutine bcast_send4D_23


subroutine bcast_recv1D_3(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 1D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx
integer :: iid


lx=size(param,1)-4

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx+4), &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv1D_3


subroutine bcast_recv1D_23_2(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 1D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx
integer :: iid


lx=size(param,1)-4

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx+4), &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv1D_23_2


subroutine bcast_recv1D_23_3(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 1D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx
integer :: iid


lx=size(param,1)-4

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx+4), &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv1D_23_3


subroutine bcast_recv2D_3(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx2,lx3
integer :: iid


lx2=size(paramtrim,1)
lx3=size(paramtrim,2)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx2*lx3, &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv2D_3


subroutine bcast_recv2D_23(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx2,lx3
integer :: iid


lx2=size(paramtrim,1)
lx3=size(paramtrim,2)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx2*lx3, &
  mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv2D_23


subroutine bcast_recv3D_3(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!> note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)    
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*lx2*lx3, &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_3


subroutine bcast_recv3D_23(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!> note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)   
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*lx2*lx3, &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_23


subroutine bcast_recv3D_x3i_3(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

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

end subroutine bcast_recv3D_x3i_3


subroutine bcast_recv3D_x3i_23(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!>note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)
lx2=size(paramtrim,2)
lx3=size(paramtrim,3)-1  ! `lx3` is an interfaced quantity

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*lx2*(lx3+1), &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_x3i_23


subroutine bcast_recv3D_x2i_23(paramtrim,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(:,:,:), intent(out) :: paramtrim
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!>note here that paramtrim does not have ghost cells
lx1=size(paramtrim,1)
lx2=size(paramtrim,2)-1    !x2 is the interfaced direction here
lx3=size(paramtrim,3)

!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(paramtrim,lx1*(lx2+1)*lx3, &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_x2i_23


subroutine bcast_recv3D_ghost_3(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!> note here that param has ghost cells
lx1=size(param,1)-4    
lx2=size(param,2)-4
lx3=size(param,3)-4


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_ghost_3


subroutine bcast_recv3D_ghost_23(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!!
!! SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!!
!! THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
!! GHOST CELLS!

real(wp), dimension(-1:,-1:,-1:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3
integer :: iid

!> note here that param has ghost cells
lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!> WORKERS RECEIVE THE IC DATA FROM ROOT
call mpi_recv(param,(lx1+4)*(lx2+4)*(lx3+4), &
               mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

end subroutine bcast_recv3D_ghost_23


subroutine bcast_recv4D_3(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,isp

lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!WORKERS RECEIVE THE IC DATA FROM ROOT
do isp=1,lsp
  call mpi_recv(param(:,:,:,isp),(lx1+4)*(lx2+4)*(lx3+4), &
                 mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end do

end subroutine bcast_recv4D_3


subroutine bcast_recv4D_23(param,tag)
!! THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
!! GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
!-------
!-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
!-------
!-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
!-------GHOST CELLS!
!------------------------------------------------------------

real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: param
integer, intent(in) :: tag

integer :: lx1,lx2,lx3,isp
real(wp), dimension(-1:size(param,1)-2,1:size(param,2)-4,1:size(param,3)-4) :: paramtmp

lx1=size(param,1)-4
lx2=size(param,2)-4
lx3=size(param,3)-4


!WORKERS RECEIVE THE IC DATA FROM ROOT
do isp=1,lsp
  call mpi_recv(paramtmp,(lx1+4)*lx2*lx3, &
                 mpi_realprec,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  param(-1:lx1+2,1:lx2,1:lx3,isp)=paramtmp(-1:lx1+2,1:lx2,1:lx3)
end do

end subroutine bcast_recv4D_23

end module mpimod
