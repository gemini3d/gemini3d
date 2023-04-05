submodule (mpimod) mpihalo

use mpi_f08, only: mpi_waitall, mpi_isend, mpi_irecv, MPI_REQUEST, MPI_STATUS, MPI_PROC_NULL

implicit none (type, external)

contains

!! HALOS AN ARRAY THAT IS SPLIT ALONG THE 2 AND 3 RANKS
!subroutine halo_23(param,lhalo,tag,isperiodic)
!real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
!integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (probably 1 or 2)
!integer, intent(in) :: tag
!logical, intent(in) :: isperiodic
module procedure halo_23
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

  integer :: lx1,lx2,lx3
  integer :: idleft,idright,idup,iddown
  integer :: i2,i3
  type(MPI_REQUEST) :: requests(4)
  type(MPI_STATUS) :: statuses(4)
  real(wp), allocatable, dimension(:,:,:) :: buffer31,buffer32,buffer33,buffer34
  real(wp), allocatable, dimension(:,:,:) :: buffer21,buffer22,buffer23,buffer24
  logical :: x2begin,x3begin,x2end,x3end
  !! to store info about whether we are the first or last process in a direction

  !COMPUTE SIZES, JUST IN CASE.
  lx1=size(param,1)-4
  lx2=size(param,2)-4
  lx3=size(param,3)-4

  !IDENTIFY MY NEIGHBORS IN X3
  x3begin=.false.
  x3end=.false.

  i3=mpi_cfg%myid3-1
  i2=mpi_cfg%myid2
  if (i3==-1) then
    !! global boundary to my left, assume periodic
    i3=mpi_cfg%lid3-1
    !! lid3-1 is the last process in x3 on the process grid
    x3begin=.true.
  end if
  idleft=grid2ID(i2,i3)
  if (x3begin .and. .not.(isperiodic)) then
    !! we are flagged as not wanting periodic boundaries so do nothing (overwrite idleft to send to NULL process
    idleft=MPI_PROC_NULL
  end if

  i3=mpi_cfg%myid3+1
  i2=mpi_cfg%myid2
  if (i3==mpi_cfg%lid3) then
    !! global boundary to my right, assume periodic
    i3=0
    x3end=.true.
  end if
  idright=grid2ID(i2,i3)
  !! convert the location on process grid into a flat processed ID, The process grid is
  !! visualized as lid2,lid3 in terms of index order (e.g. the i2 index cycles more quickly
  if (x3end .and. .not.(isperiodic)) then
    idright=MPI_PROC_NULL
  end if

  !IDENTIFY MY NEIGHBORING PROCESSES IN X2
  x2begin=.false.
  x2end=.false.

  i3=mpi_cfg%myid3
  i2=mpi_cfg%myid2-1
  if (i2==-1) then
    !! global boundary downward, assume periodic
    i2=mpi_cfg%lid2-1
    x2begin=.true.
  end if
  iddown=grid2ID(i2,i3)
  if (x2begin) then
    !! never assume periodic in the x2-direction
    iddown=MPI_PROC_NULL
  end if

  i3=mpi_cfg%myid3
  i2=mpi_cfg%myid2+1
  if (i2==mpi_cfg%lid2) then
    !! global boundary upward, assume periodic
    i2=0
    x2end=.true.
  end if
  idup=grid2ID(i2,i3)
  !! convert to process ID
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
  if (.not. (x3begin .and. x3end)) then
    !! make sure we actually need to pass in this direction, viz. we aren't both the beginning and thend
    buffer31=param(-1:lx1+2,1:lx2,1:lhalo)     !x1 ghost cells to be included
    call mpi_isend(buffer31,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft,tag,MPI_COMM_WORLD, requests(1))

    buffer32=param(-1:lx1+2,1:lx2,lx3+1-lhalo:lx3)
    call mpi_isend(buffer32,(lx1+4)*(lx2)*lhalo,mpi_realprec, &
                      idright,tag,MPI_COMM_WORLD, requests(2))

    call mpi_irecv(buffer33,(lx1+4)*(lx2)*lhalo,mpi_realprec,idright, &
                      tag,MPI_COMM_WORLD, requests(3))

    call mpi_irecv(buffer34,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft, &
                            tag,MPI_COMM_WORLD, requests(4))

    call mpi_waitall(4,requests,statuses)

    if (idright/=MPI_PROC_NULL)  then    !only overwrite the cells if we didn't do a null receive
      param(-1:lx1+2,1:lx2,lx3+1:lx3+lhalo)=buffer33    !can't copy out buffers until we know the messages have been received
    end if
    if (idleft/=MPI_PROC_NULL) then
      param(-1:lx1+2,1:lx2,1-lhalo:0)=buffer34
    end if
  else if (isperiodic) then   ! there is the possibility on a periodic grid with lid3=1 that we still need to enforce periodic conditions.  Note that this condition is separate from whether or not the workers is the first and last x3 worker...
    !param(-1:lx1+2,1:lx2,1-lhalo:0)=param(-1:lx1+2,1:lx2,lx3-1:lx3)
    !param(-1:lx1+2,1:lx2,lx3+1:lx3+lhalo)=param(-1:lx1+2,1:lx3,1:2)  !yikes, wrong size+seg fault if lhalo=1???
    param(-1:lx1+2,1:lx2,1-lhalo:0)=param(-1:lx1+2,1:lx2,lx3-(lhalo-1):lx3)
    param(-1:lx1+2,1:lx2,lx3+1:lx3+lhalo)=param(-1:lx1+2,1:lx2,1:lhalo)
  end if

  !EXCHANGE MESSAGES IN THE X2 DIRECTION OF THE PROCESS GRID
  if (.not. (x2begin .and. x2end)) then
    buffer21=param(-1:lx1+2,1:lhalo,1:lx3)
    call mpi_isend(buffer21,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown,tag,MPI_COMM_WORLD, requests(1))

    buffer22=param(-1:lx1+2,lx2+1-lhalo:lx2,1:lx3)
    call mpi_isend(buffer22,(lx1+4)*(lx3)*lhalo,mpi_realprec, &
                      idup,tag,MPI_COMM_WORLD, requests(2))

    call mpi_irecv(buffer23,(lx1+4)*(lx3)*lhalo,mpi_realprec,idup,&
                      tag,MPI_COMM_WORLD, requests(3))

    call mpi_irecv(buffer24,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown, &
                            tag,MPI_COMM_WORLD, requests(4))

    call mpi_waitall(4,requests,statuses)

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
end procedure halo_23


!> This halos a 4D array (4-dimension is species dimension for GEMINI)
!subroutine halo_allspec_23(param,lhalo,tag,isperiodic)
!real(wp), dimension(-1:,-1:,-1:,:), intent(inout) :: param
!integer, intent(in) :: lhalo    !number of surrounding grid points to halo with (probably 1 or 2)
!integer, intent(in) :: tag
!logical, intent(in) :: isperiodic
module procedure halo_allspec_23
  integer :: isp,lsp

  lsp=size(param,4)

  do isp=1,lsp
    call halo(param(:,:,:,isp),lhalo,tag,isperiodic)
  end do
end procedure halo_allspec_23


!subroutine halo_end_23(param,paramend,paramtop,paramcorner,tag)
!real(wp), dimension(:,:,:), intent(inout) :: param
!real(wp), dimension(:,:), intent(inout) :: paramend
!real(wp), dimension(:,:), intent(inout) :: paramtop
!real(wp), dimension(:), intent(inout) :: paramcorner
!integer, intent(in) :: tag
module procedure halo_end_23
  !! GENERIC HALOING ROUTINE WHICH PASSES THE BEGINNING OF THE
  !! SLAB TO ITS LEFTWARD (IN X3) NEIGHBOR SO THAT X3 INTEGRATIONS
  !! CAN BE DONE PROPERLY.  PRESENTLY THIS IS JUST USED IN MAGCALC

  integer :: lx1,lx2,lx3
  integer :: idleft,idright,iddown,idup,iddownleft,idupright
  integer :: i2,i3

  type(MPI_REQUEST) :: requests(2)
  type(MPI_STATUS) :: statuses(4)

  logical :: x2begin,x2end,x3begin,x3end,downleft,upright
  real(wp), dimension(:,:), allocatable :: buffer
  real(wp), dimension(:), allocatable :: buf_corner

  !> system sizes based off of input data
  lx1=size(param,1)
  lx2=size(param,2)
  lx3=size(param,3)

  !> identify neighbors in x3, we send our data "left" (i3-1) and receive from our "right" (i3+1)
  x3begin=.false.
  x3end=.false.

  i3=mpi_cfg%myid3-1
  i2=mpi_cfg%myid2
  if (i3==-1) then
    !! global boundary to my left, assume periodic
    i3=mpi_cfg%lid3-1
    !! lid3-1 is the last process in x3 on the process grid
    x3begin=.true.
  end if
  idleft=grid2ID(i2,i3)
  if (x3begin) then
    !! we are flagged as not wanting periodic boundaries so do nothing (overwrite idleft to send to NULL process
    idleft=MPI_PROC_NULL
  end if

  i3=mpi_cfg%myid3+1
  i2=mpi_cfg%myid2
  if (i3==mpi_cfg%lid3) then
    !! global boundary to my right, assume periodic
    i3=0
    x3end=.true.
  end if
  idright=grid2ID(i2,i3)
  !! convert the location on process grid into a flat processed ID, The process grid is
  !! visualized as lid2,lid3 in terms of index order (e.g. the i2 index cycles more quickly
  if (x3end) then
    idright=MPI_PROC_NULL
  end if

  !> identify x2 neighbor processes
  x2begin=.false.
  x2end=.false.

  i3=mpi_cfg%myid3
  i2=mpi_cfg%myid2-1
  if (i2==-1) then
    !! global boundary downward, assume periodic
    i2=mpi_cfg%lid2-1
    x2begin=.true.
  end if
  iddown=grid2ID(i2,i3)
  if (x2begin) then
    !! never assume periodic in the x2-direction
    iddown=MPI_PROC_NULL
  end if

  i3=mpi_cfg%myid3
  i2=mpi_cfg%myid2+1
  if (i2==mpi_cfg%lid2) then
    !! global boundary upward, assume periodic
    i2=0
    x2end=.true.
  end if
  idup = grid2ID(i2,i3)
  !! convert to process ID
  if (x2end) then
    idup=MPI_PROC_NULL
  end if

  !> need to identify "corner neighbor" processes
  downleft=.false.
  upright=.false.

  if (.not. (x2begin .or. x3begin)) then
    !! no down/left corner point if we reside on a min x2,3 edge
    i3=mpi_cfg%myid3-1
    i2=mpi_cfg%myid2-1
    iddownleft=grid2ID(i2,i3)
    downleft=.true.
  else
    iddownleft=MPI_PROC_NULL
  end if
  if (.not. (x2end .or. x3end)) then
    i3=mpi_cfg%myid3+1
    i2=mpi_cfg%myid2+1
    idupright=grid2ID(i2,i3)
    upright=.true.
  else
    idupright=MPI_PROC_NULL
  end if

  !> data passing in x3, if appropriate
  if (.not. (x3begin .and. x3end)) then
    !! for singleton process grid along x3; do not want to try to send to self...
    !! make sure we actually need to pass in this direction, viz. we aren't both the beginning and thend

    !> force contiguous in memory
    allocate(buffer(lx1,lx2))
    buffer=param(:,:,1)
    call mpi_isend(buffer,lx1*lx2,mpi_realprec,idleft,tag,MPI_COMM_WORLD, requests(1))

    call mpi_irecv(buffer,lx1*lx2,mpi_realprec,idright,tag,MPI_COMM_WORLD, requests(2))
    paramend = buffer

    call mpi_waitall(2,requests,statuses)
    deallocate(buffer)
  end if

  !> data passing in x2, if appropriate
  if (.not. (x2begin .and. x2end)) then
    !! for singleton process grid along x2; dont' send to self

    !> force contiguous in memory
    allocate(buffer(lx1,lx3))
    buffer = param(:,1,:)
    call mpi_isend(buffer,lx1*lx3,mpi_realprec,iddown,tag,MPI_COMM_WORLD, requests(1))

    call mpi_irecv(buffer,lx1*lx3,mpi_realprec,idup,tag,MPI_COMM_WORLD, requests(2))
    paramtop = buffer

    call mpi_waitall(2,requests,statuses)
    deallocate(buffer)
  end if

  !> corner data passing if necessary
  if (.not. (x2begin .and. x2end .and. x3begin .and. x3end)) then
    !! single process "corner" case, lol; don't send to self

    !> force data into a contiguous buffer
    allocate(buf_corner(lx1))
    buf_corner=param(:,1,1)
    call mpi_isend(buf_corner,lx1,mpi_realprec,iddownleft,tag,MPI_COMM_WORLD, requests(1))

    call mpi_irecv(buf_corner,lx1,mpi_realprec,idupright,tag,MPI_COMM_WORLD, requests(2))
    paramcorner = buf_corner

    call mpi_waitall(2,requests,statuses)
    deallocate(buf_corner)
  end if

  !zero out ghost cells if past the end of the full simulation grid
  if (mpi_cfg%myid2==mpi_cfg%lid2-1) paramtop=0
  !! add nothing on the end since no one is passing leftward to me, FIXME: need to account for periodic???
  if (mpi_cfg%myid3==mpi_cfg%lid3-1) paramend=0
  !! zero out the data at the end of the grid
  if (mpi_cfg%myid2==mpi_cfg%lid2-1 .or. mpi_cfg%myid3==mpi_cfg%lid3-1) paramcorner=0
end procedure halo_end_23


end submodule mpihalo
