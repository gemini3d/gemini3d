module mpimod


!NOTES:
! - Need to consider overloading routines as send_ghost and send_noghost so that
! it is more clear what the structure of the input arrays should be. 

use phys_consts, only : lsp   !code needs to know how many species are being used.
use mpi, only: mpi_init, mpi_comm_world, mpi_double_precision, mpi_status_size, mpi_status_ignore

implicit none

private :: lsp
private :: mpi_init, mpi_comm_world, mpi_double_precision, mpi_status_size, mpi_status_ignore

!NOW A LIST OF TAGS SO THESE DO NOT NEED TO BE EMBEDDED IN EACH SUBROUTINE
integer, parameter :: tagns=2,tagvs1=3,tagTs=4    !root/workers input routines.  also output routines for root/worker
integer, parameter :: tagJ1=6,tagJ2=7,tagJ3=8,tagv2=9,tagv3=10    !output root/worker routines, main program, potential_comm
integer, parameter :: tagvs3BC=100,tagnsBC=101,tagrhovs1BC=102,tagrhoesBC=103    !used in the advection boundary conditions
integer, parameter :: tagvs1BC=1000,tagvs2BC=1001!,tagvs3BC=1002    !used in the compression solution

!THESE MESSAGES ARE USED IN ELECTRODYNAMICS MODULE
integer, parameter :: tagE1=11,tagE2=12,tagE3=13
integer, parameter :: tagsigP=16,tagsigH=17,tagsig0=18,tagincap=19,tagv2pol=20,tagv3pol=21
integer, parameter :: tagDE2Dt=22,tagDE3Dt=23,tagflagdirich=24    !unused in present version
integer, parameter :: tagvn2=25,tagvn3=26,tagB1=27                !for passing/gathering full-grid winds

!IN THE MAIN PROGRAM
integer, parameter :: tagx3=1,tagdt=5
integer, parameter :: tagx1=27,tagx2=28

!IN THE GRID MODULE
integer, parameter :: tagh1=29,tagh2=30,tagh3=31
integer, parameter :: tagglat=32,tagglon=33,tagalt=34
integer, parameter :: taglx1=35,taglx2=36,taglx3=37,taglx3all=38
integer, parameter :: tagBmag=39,taginc=40,tagnull=41
integer, parameter :: tageunit1=42,tageunit2=43,tageunit3=44,tager=45,tagetheta=46,tagephi=47
integer, parameter :: tagr=56,tagtheta=57,tagphi=58

!IN THE NEUTRAL MODULE
integer, parameter :: taglrho=48,taglz=49
integer, parameter :: tagdnO=50,tagdnN2=51,tagdnO2=52,tagdTn=53,tagdvnrho=54,tagdvnz=55
integer, parameter :: tagly=69

!FOR DEALING WITH PRECIPITATION BOUNDARY CONDITIONS MODULE
integer, parameter :: tagllat=59,tagllon=60,tagmlat=61,tagmlon=62,tagQp=63,tagE0p=64

!FOR DEALING WITH THE ELECTRIC FIELD BOUNDARY CONDITIONS
integer, parameter :: tagE0xp=65,tagE0yp=66,tagE0xi=67,tagE0yi=68

!FOR DISTRIBUTING PART OF THE ELECTRODYNAMICS CALCULATIONS
integer, parameter :: tagsrc=69,tagSigPint2=70,tagSigPint3=71,tagSigHint=72,tagincapint=73,tagv2electro=74,tagv3electro=75
integer, parameter :: tagE01=76,tagE02=77,tagE03=78,tagVminx1=79,tagVmaxx1=80

!THESE ARE USED IN MAGCALC.F90 PROGRAM
integer, parameter :: tagBr=81,tagBtheta=82, tagBphi=83

!VARIABLES REUSED BY ALL WORKERS AND USING MODULES
integer, protected :: myid,lid    !no external procedure should mess with these (but they need to be able to read them)
integer :: ierr                   !using procedures need to be able to overwrite this to prevent seg. faults (or something)


!VARIABLES RELATED TO PROCESS GRID (IF USED)
integer, protected :: lid2,lid3,myid2,myid3


!THESE INTERFACES OVERLOAD THE MPI GATHER,BROADCAST SUBROUTINES FOR ARRAYS OF DIFFERENT RANKS.
interface gather_recv
  module procedure gather_recv2D, gather_recv3D, gather_recv4D
end interface gather_recv

interface gather_send
  module procedure gather_send2D, gather_send3D, gather_send4D
end interface gather_send

interface bcast_send
  module procedure bcast_send1D, bcast_send2D, bcast_send3D, bcast_send4D
end interface bcast_send

interface bcast_recv
  module procedure bcast_recv1D, bcast_recv2D, bcast_recv3D, bcast_recv4D
end interface bcast_recv


contains


  subroutine mpisetup()

    !------------------------------------------------------------
    !-------THIS SUBROUTINE INITIALIZES MODULE MPI VARIABLES FOR A 
    !-------WORKER.  THIS CURRENTLY IS UNUSED AS IT HAS NOT BEEN
    !------FULLY IMPLEMENTED IN THIS VERSINO OF THE CODE.
    !------------------------------------------------------------ 

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,lid,ierr)

    !INITIALIZE WITH ONLY PARALLELIZING IN X3, GRIDDING FUNCTION MAY CHANGE
    !THIS, IF CALLED.
    lid2=1
    lid3=lid

  end subroutine mpisetup


  subroutine mpigrid(lx2all,lx3all)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE DEFINES A PROCESS GRID, IF REQUIRED 
    !-------IT IS CURRENTLY NOT USED BUT KEPT HERE FOR FUTURE DEVELOPMENT
    !------------------------------------------------------------ 

    integer, intent(in) :: lx2all,lx3all

    lid2=1
    lid3=lid
    do while( ((lid3/2)*2==lid3) .and. (lid3-lid2>lid3 .or. lid3-lid2>lid2) .and. &     !also must insure that lx3 is divisible by lid3 and lx2 by lid2
             lx3all/(lid3/2)*(lid3/2)==lx3all .and. lx2all/(lid2*2)*(lid2*2)==lx2all)
      lid3=lid3/2
      lid2=lid2*2
    end do


    !THIS PROCESS' LOCATION ON THE GRID
    myid3=myid/lid2
    myid2=myid-myid3*lid2

    write(*,*) 'Proposed process grid is x2 by x3 size (in number of processes):  ',lid2,' by ',lid3
    write(*,*) 'Process:  ',myid,' is at location:  ',myid2,myid3,' on the process grid'

  end subroutine mpigrid


  subroutine mpibreakdown()

    !------------------------------------------------------------
    !-------THIS SUBROUTINE SHUTS DOWN MPI
    !------------------------------------------------------------

    call mpi_finalize(ierr)

  end subroutine mpibreakdown


  subroutine halo(param,lhalo,tag)

    !------------------------------------------------------------
    !-------GENERIC HALOING ROUTINE FOR FILLING GHOST CELLS.  CAN
    !-------BE USED TO SET BOUNDARY CONDITIONS OR PREPARE ARRAYS
    !-------FOR FINITE DIFFERENCING, ETC.  OBVIOUSLY ARRAYS INCLUDE
    !-------GHOST CELLS.  ARRAYS SHOULD HAVE SPECIES DIMENSION
    !-------REMOVED BEFORE PASSAGE INTO THIS SUBROUTINE.  THIS
    !-------ROUTINE WORKS FOR A 3D ARRAY, AND BY DEFAULT ENFORCES
    !-------PERIODIC BOUNDARY CONDITIONS.  IF APERIODIC CONDITIONS
    !-------ARE NEEDED YOU MUST OVERWRITE THE GLOBAL BOUNDARIES AFTER
    !-------YOU CALL HALO.
    !-------
    !-------THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:), intent(inout) :: param
    integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (1 or 2 only)
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3,ihalo
    integer :: idleft,idright

    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: statuses
    integer :: tmpreq

    real(8) :: tstart,tfin

    lx1=size(param,1)-4
    lx2=size(param,2)-4
    lx3=size(param,3)-4


    !IDENTIFY MY NEIGHBORS
    idleft=myid-1
    idright=myid+1


    !SCREEN FOR GLOBAL BOUNDARIES, ASSUME PERIODIC (MUST BE OVERWRITTEN LATER IF
    !YOU ARE USING ANOTHER TYPE OF BOUNDARY
    if (idleft==-1) then
      idleft=lid-1
!      idleft=MPI_PROC_NULL    !if you wanted to default to aperiodic you could do this...
    end if
    if (idright==lid) then
      idright=0
!      idright=MPI_PROC_NULL
    end if


    call mpi_isend(param(:,:,1:lhalo),(lx1+4)*(lx2+4)*lhalo,MPI_DOUBLE_PRECISION,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(1)=tmpreq
    call mpi_isend(param(:,:,lx3+1-lhalo:lx3),(lx1+4)*(lx2+4)*lhalo,MPI_DOUBLE_PRECISION, &
                      idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(2)=tmpreq
    call mpi_irecv(param(:,:,lx3+1:lx3+lhalo),(lx1+4)*(lx2+4)*lhalo,MPI_DOUBLE_PRECISION,idright, &
                      tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(3)=tmpreq
    call mpi_irecv(param(:,:,1-lhalo:0),(lx1+4)*(lx2+4)*lhalo,MPI_DOUBLE_PRECISION,idleft, &
                            tag,MPI_COMM_WORLD,tmpreq,ierr)
    requests(4)=tmpreq

    call mpi_waitall(4,requests,statuses,ierr)

  end subroutine halo


!   LEAVE THIS IN FOR FUTURE DEV. - HALOS AN ARRAY THAT IS SPLIT ALONG THE 2 AND
!   3 RANKS
!
!  subroutine halo23(param,lhalo,tag)
!
!    !------------------------------------------------------------
!    !-------GENERIC HALOING ROUTINE FOR FILLING GHOST CELLS.  CAN
!    !-------BE USED TO SET BOUNDARY CONDITIONS OR PREPARE ARRAYS
!    !-------FOR FINITE DIFFERENCING, ETC.  OBVIOUSLY ARRAYS INCLUDE
!    !-------GHOST CELLS.  ARRAYS SHOULD HAVE SPECIES DIMENSION
!    !-------REMOVED BEFORE PASSAGE INTO THIS SUBROUTINE.  THIS
!    !-------ROUTINE WORKS FOR A 3D ARRAY.
!    !-------
!    !-------THIS SUBROUTINE SHOULD BE CONSIDERED DISTINCT IN FUNCTIONALITY
!    !-------FROM THE SUBROUTINE WHICH COMPUTES BOUNDARY CONDITIONS FOR
!    !-------THE GLOBAL GRID.  I.E., IF CALLED ON A WORKER THAT ABUTS THE
!    !-------GLOBAL BOUNDARY IT WILL DO *NOTHING*.
!    !-------
!    !-------THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS
!    !-------
!    !-------THIS VERSION ALSO ASSUMES A PROCESS GRID HAS BEEN DEFINED
!    !-------AND THAT PASSING NEEDS TO BE DONE IN X2 AND X3
!    !------------------------------------------------------------
!
!    real(8), dimension(-1:,-1:,-1:), intent(inout) :: param
!    integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (probably 1 or 2)
!    integer, intent(in) :: tag
!
!    integer :: lx1,lx2,lx3,ihalo
!    integer :: idleft,idright,idup,iddown
!    integer :: i2,i3
!
!    integer, dimension(4) :: requests
!    integer, dimension(MPI_STATUS_SIZE,4) :: statuses
!    integer :: tmpreq
!
!    real(8) :: tstart,tfin
!
!    lx1=size(param,1)-4
!    lx2=size(param,2)-4
!    lx3=size(param,3)-4
!
!
!    !IDENTIFY MY NEIGHBORS IN X3
!    i3=myid3-1
!    if (i3==-1) then    !global boundary to my left
!      i3=lid3
!    end if
!    i2=myid2
!    idleft=(i3-1)*lid2+i2
!    i3=myid3+1
!    if (i3==lid3) then    !global boundary to my right
!      i3=1
!    end if
!    i2=myid2
!    idright=(i3-1)*lid2+i2
!
!
!    !IDENTIFY MY NEIGHBORING PROCESSES IN X2
!    i3=myid3
!    i2=myid2-1
!    if (i2==-1) then    !global boundary downward
!      i2=lid2
!    end if
!    iddown=(i3-1)*lid2+i2
!    i3=myid3
!    i2=myid2+1
!    if (i2==lid2) then     !global boundary upward
!      i2=1
!    end if
!    idup=(i3-1)*lid2+i2
!
!
!    !ALLOCATE SPACE TO BUFFER MESSAGES (SINCE USING ASYNCHRONOUS MPI COMMANDS)
!    allocate(buffer31(-1:lx1+2,1:lx2,lhalo),buffer32(-1:lx1+2,1:lx2,lhalo),buffer33(-1:lx1+2,1:lx2,lhalo),buffer34(-1:lx1+2,1:lx2,lhalo))
!    allocate(buffer21(-1:lx1+2,1:lx2,lhalo),buffer22(-1:lx1+2,1:lx2,lhalo),buffer23(-1:lx1+2,1:lx2,lhalo),buffer24(-1:lx1+2,1:lx2,lhalo))
!
!
!    !EXCHANGE MESSAGES IN THE X3-DIRECTION OF THE PROCESS GRID
!    buffer31=param(-1:lx1+2,1:lx2,1:lhalo)     !x1 ghost cells to be included
!    call mpi_isend(buffer31,(lx1+4)*(lx2)*lhalo,MPI_DOUBLE_PRECISION,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(1)=tmpreq
!
!    buffer32=param(-1:lx1+2,1:lx2,lx3+1-lhalo:lx3)
!    call mpi_isend(buffer32,(lx1+4)*(lx2)*lhalo,MPI_DOUBLE_PRECISION, &
!                      idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(2)=tmpreq
!
!    call mpi_irecv(buffer33,(lx1+4)*(lx2)*lhalo,MPI_DOUBLE_PRECISION,idright, &
!                      tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(3)=tmpreq
!
!    call mpi_irecv(buffer34,(lx1+4)*(lx2)*lhalo,MPI_DOUBLE_PRECISION,idleft, &
!                            tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(4)=tmpreq
!
!    call mpi_waitall(4,requests,statuses,ierr)
!    param(-1:lx1+2,1:lx2,lx3+1:lx3+lhalo)=buffer33    !can't copy out buffers until we know the messages have been received
!    param(-1:lx1+2,1:lx2,1-lhalo:0)=buffer34
!
!
!    !EXCHANGE MESSAGES IN THE X2 DIRECTION OF THE PROCESS GRID
!    buffer21=param(-1:lx1+2,1:lhalo,1:lx3)
!    call mpi_isend(buffer21,(lx1+4)*(lx3)*lhalo,MPI_DOUBLE_PRECISION,iddown,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(1)=tmpreq
!
!    buffer22=param(-1:lx1+2,lx2+1-lhalo:lx2,1:lx3)
!    call mpi_isend(buffer22,(lx1+4)*(lx3)*lhalo,MPI_DOUBLE_PRECISION, &
!                      idup,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(2)=tmpreq
!
!    call mpi_irecv(buffer23,(lx1+4)*(lx3)*lhalo,MPI_DOUBLE_PRECISION,idup,&
!                      tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(3)=tmpreq
!
!    call mpi_irecv(buffer24,(lx1+4)*(lx3)*lhalo,MPI_DOUBLE_PRECISION,iddown, &
!                            tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(4)=tmpreq
!    
!    call mpi_waitall(4,requests,statuses,ierr)
!    param(-1:lx1+2,lx2+1:lx2+lhalo,1:lx3)=buffer23    !clear to copy out buffers
!    param(-1,lx1+2,1-lhalo:0,1:lx3)=buffer24
!
!
!    !CLEAR OUT BUFFER VARIABLES
!    deallocate(buffer31,buffer32,buffer33,buffer34)
!    deallocate(buffer21,buffer22,buffer23,buffer24)
!
!  end subroutine halo23


  subroutine gather_recv2D(paramtrim,tag,paramtrimall)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
    !-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
    !-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
    !-------
    !-------THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
    !-------
    !-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
    !-------ANY GHOST CELLS!!!!
    !------------------------------------------------------------

    real(8), dimension(:,:), intent(in) :: paramtrim
    integer, intent(in) :: tag
    real(8), dimension(:,:), intent(out) :: paramtrimall

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
                    MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end do

  end subroutine gather_recv2D


  subroutine gather_recv3D(paramtrim,tag,paramtrimall)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE GATHERS DATA FROM ALL WORKERS ONTO
    !-------A FULL-GRID ARRAY ON THE ROOT PROCESS (PRESUMABLY FOR
    !-------OUTPUT OR SOME ELECTRODYNAMIC CALCULATION, PERHAPS.
    !-------
    !-------THIS SUBROUTINE IS TO BE CALLED BY ROOT TO DO GATHER
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------ANY GHOST CELLS!!!!
    !------------------------------------------------------------

    real(8), dimension(:,:,:), intent(in) :: paramtrim
    integer, intent(in) :: tag
    real(8), dimension(:,:,:), intent(out) :: paramtrimall

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
                    MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end do

  end subroutine gather_recv3D


  subroutine gather_recv4D(param,tag,paramall)

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

    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: param
    integer, intent(in) :: tag
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: paramall

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
                      MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      end do
    end do

  end subroutine gather_recv4D


  subroutine gather_send2D(paramtrim,tag)

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

    real(8), dimension(:,:), intent(in) :: paramtrim
    integer, intent(in) :: tag

    integer :: lx2,lx3


    lx2=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
    lx3=size(paramtrim,2)

    call mpi_send(paramtrim,lx2*lx3,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)

  end subroutine gather_send2D


  subroutine gather_send3D(paramtrim,tag)

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

    real(8), dimension(:,:,:), intent(in) :: paramtrim
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3


    lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
    lx2=size(paramtrim,2)
    lx3=size(paramtrim,3)

    call mpi_send(paramtrim,lx1*lx2*lx3,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)

  end subroutine gather_send3D


  subroutine gather_send4D(param,tag)

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

    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: param
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3,isp


    lx1=size(param,1)-4
    lx2=size(param,2)-4
    lx3=size(param,3)-4

    do isp=1,lsp
      call mpi_send(param(:,:,1:lx3,isp),(lx1+4)*(lx2+4)*lx3,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,ierr)
    end do

  end subroutine gather_send4D


  subroutine bcast_send1D(paramall,tag,param)

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

    real(8), dimension(-1:), intent(in) :: paramall
    integer, intent(in) :: tag
    real(8), dimension(-1:), intent(out) :: param

    integer :: lx,lxall     !local sizes
    integer :: iid,islstart,islfin


    lxall=size(paramall,1)-4
    lx=size(param,1)-4


    do iid=1,lid-1
      islstart=iid*lx+1
      islfin=islstart+lx-1

      call mpi_send(paramall(islstart-2:islfin+2),(lx+4), &
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)
    end do
    param=paramall(-1:lx+2)

  end subroutine bcast_send1D


  subroutine bcast_send2D(paramtrimall,tag,paramtrim)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
    !-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 2D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(:,:), intent(in) :: paramtrimall
    integer, intent(in) :: tag
    real(8), dimension(:,:), intent(out) :: paramtrim

    integer :: lx2,lx3
    integer :: iid,islstart,islfin


    lx2=size(paramtrim,1)    !assume this is an array which has been 'flattened' along the 1-dimension
    lx3=size(paramtrim,2)


    !ROOT BROADCASTS IC DATA TO WORKERS
    do iid=1,lid-1
      islstart=iid*lx3+1
      islfin=islstart+lx3-1

      call mpi_send(paramtrimall(:,islstart:islfin),lx2*lx3, &
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)
    end do


    !ROOT TAKES A SLAB OF DATA
    paramtrim=paramtrimall(:,1:lx3)

  end subroutine bcast_send2D


  subroutine bcast_send3D(paramtrimall,tag,paramtrim)

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

    real(8), dimension(:,:,:), intent(in) :: paramtrimall
    integer, intent(in) :: tag
    real(8), dimension(:,:,:), intent(out) :: paramtrim

    integer :: lx1,lx2,lx3
    integer :: iid,islstart,islfin


    lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
    lx2=size(paramtrim,2)
    lx3=size(paramtrim,3)


    !ROOT BROADCASTS IC DATA TO WORKERS
    do iid=1,lid-1
      islstart=iid*lx3+1
      islfin=islstart+lx3-1

      call mpi_send(paramtrimall(:,:,islstart:islfin),lx1*lx2*lx3, &
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)
    end do


    !ROOT TAKES A SLAB OF DATA
    paramtrim=paramtrimall(:,:,1:lx3)

  end subroutine bcast_send3D


  subroutine bcast_send3D_x3i(paramtrimall,tag,paramtrim)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
    !-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS, BUT ARE X3 INTERFACE QUANITTIES
    !------------------------------------------------------------

    real(8), dimension(:,:,:), intent(in) :: paramtrimall
    integer, intent(in) :: tag
    real(8), dimension(:,:,:), intent(out) :: paramtrim

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
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)     !note the +1 since thes are interfact quantities (and need to overlap b/t workers)
    end do


    !ROOT TAKES A SLAB OF DATA
    paramtrim=paramtrimall(:,:,1:lx3+1)

  end subroutine bcast_send3D_x3i


  subroutine bcast_send3D_ghost(paramall,tag,param)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
    !-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH INCLUDE GHOST CELLS
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:), intent(in) :: paramall
    integer, intent(in) :: tag
    real(8), dimension(-1:,-1:,-1:), intent(out) :: param

    integer :: lx1,lx2,lx3
    integer :: iid,islstart,islfin


    lx1=size(param,1)-4    !note here that param has ghost cells
    lx2=size(param,2)-4
    lx3=size(param,3)-4


    !ROOT BROADCASTS IC DATA TO WORKERS
    do iid=1,lid-1
      islstart=iid*lx3+1
      islfin=islstart+lx3-1

      call mpi_send(paramall(:,:,islstart-2:islfin+2),(lx1+4)*(lx2+4)*(lx3+4), &
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)
    end do


    !ROOT TAKES A SLAB OF DATA
    param=paramall(:,:,-1:lx3+2)

  end subroutine bcast_send3D_ghost


  subroutine bcast_send4D(paramall,tag,param)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE BROADCASTS DATA FROM A FULL GRID ARRAY
    !-------ON ROOT PROCESS TO ALL WORKERS' SUB-GRID ARRAYS.
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY ROOT TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: paramall
    integer, intent(in) :: tag
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: param

    integer :: lx1,lx2,lx3,isp
    integer :: iid,islstart,islfin


    lx1=size(param,1)-4
    lx2=size(param,2)-4
    lx3=size(param,3)-4


    !ROOT BROADCASTS IC DATA TO WORKERS
    do isp=1,lsp
      param(:,:,:,isp)=paramall(:,:,-1:lx3,isp)    !roots part of the data

      do iid=1,lid-1
        islstart=iid*lx3+1
        islfin=islstart+lx3-1

        call mpi_send(paramall(:,:,islstart-2:islfin+2,isp),(lx1+4)*(lx2+4)*(lx3+4), &
                   MPI_DOUBLE_PRECISION,iid,tag,MPI_COMM_WORLD,ierr)
      end do
    end do

  end subroutine bcast_send4D


  subroutine bcast_recv1D(param,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 1D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(-1:), intent(out) :: param
    integer, intent(in) :: tag

    integer :: lx
    integer :: iid


    lx=size(param,1)-4

    !WORKERS RECEIVE THE IC DATA FROM ROOT
    call mpi_recv(param,(lx+4), &
                   MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  end subroutine bcast_recv1D


  subroutine bcast_recv2D(paramtrim,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(:,:), intent(out) :: paramtrim
    integer, intent(in) :: tag

    integer :: lx2,lx3
    integer :: iid


    lx2=size(paramtrim,1)
    lx3=size(paramtrim,2)


    !WORKERS RECEIVE THE IC DATA FROM ROOT
    call mpi_recv(paramtrim,lx2*lx3, &
                   MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  end subroutine bcast_recv2D


  subroutine bcast_recv3D(paramtrim,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(:,:,:), intent(out) :: paramtrim
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3
    integer :: iid


    lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
    lx2=size(paramtrim,2)
    lx3=size(paramtrim,3)


    !WORKERS RECEIVE THE IC DATA FROM ROOT
    call mpi_recv(paramtrim,lx1*lx2*lx3, &
                   MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  end subroutine bcast_recv3D


  subroutine bcast_recv3D_x3i(paramtrim,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(:,:,:), intent(out) :: paramtrim
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3
    integer :: iid


    lx1=size(paramtrim,1)    !note here that paramtrim does not have ghost cells
    lx2=size(paramtrim,2)
    lx3=size(paramtrim,3)-1    !this is an x3i quantity


    !WORKERS RECEIVE THE IC DATA FROM ROOT
    call mpi_recv(paramtrim,lx1*lx2*(lx3+1), &
                   MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  end subroutine bcast_recv3D_x3i


  subroutine bcast_recv3D_ghost(param,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 3D ARRAYS WHICH DO NOT INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:), intent(out) :: param
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3
    integer :: iid


    lx1=size(param,1)-4    !note here that param has ghost cells
    lx2=size(param,2)-4
    lx3=size(param,3)-4


    !WORKERS RECEIVE THE IC DATA FROM ROOT
    call mpi_recv(param,(lx1+4)*(lx2+4)*(lx3+4), &
                   MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  end subroutine bcast_recv3D_ghost


  subroutine bcast_recv4D(param,tag)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE RECEIVES BROADCAST DATA FROM A FULL 
    !-------GRID ARRAY ON ROOT PROCESS TO WORKERS' SUB-GRID ARRAYS. 
    !-------
    !-------SUBROUTINE IS TO BE CALLED BY WORKERS TO DO A BROADCAST
    !-------
    !-------THIS VERSION WORKS ON 4D ARRAYS WHICH INCLUDE
    !-------GHOST CELLS!
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: param
    integer, intent(in) :: tag

    integer :: lx1,lx2,lx3,isp


    lx1=size(param,1)-4
    lx2=size(param,2)-4
    lx3=size(param,3)-4


    !WORKERS RECEIVE THE IC DATA FROM ROOT
    do isp=1,lsp
      call mpi_recv(param(:,:,:,isp),(lx1+4)*(lx2+4)*(lx3+4), &
                     MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end do

  end subroutine bcast_recv4D

end module mpimod
