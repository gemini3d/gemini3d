submodule (mpimod) mpihalo

implicit none

contains

module procedure halo
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

end procedure halo


module procedure halo_end
!! GENERIC HALOING ROUTINE WHICH PASSES THE BEGINNING OF THE
!! SLAB TO ITS LEFTWARD (IN X3) NEIGHBOR SO THAT X3 INTEGRATIONS
!! CAN BE DONE PROPERLY
!!
!! THIS VERSION USES ASYNC COMM WITHOUT SWITCH STATEMENTS.  IT
!! ALSO ASSUMES THAT THE ARRAYS INVOLVED DO HAVE GHOST CELLS

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

end procedure halo_end


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
!    real(wp), dimension(-1:,-1:,-1:), intent(inout) :: param
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
!    real(wp) :: tstart,tfin
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
!    call mpi_isend(buffer31,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(1)=tmpreq
!
!    buffer32=param(-1:lx1+2,1:lx2,lx3+1-lhalo:lx3)
!    call mpi_isend(buffer32,(lx1+4)*(lx2)*lhalo,mpi_realprec, &
!                      idright,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(2)=tmpreq
!
!    call mpi_irecv(buffer33,(lx1+4)*(lx2)*lhalo,mpi_realprec,idright, &
!                      tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(3)=tmpreq
!
!    call mpi_irecv(buffer34,(lx1+4)*(lx2)*lhalo,mpi_realprec,idleft, &
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
!    call mpi_isend(buffer21,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(1)=tmpreq
!
!    buffer22=param(-1:lx1+2,lx2+1-lhalo:lx2,1:lx3)
!    call mpi_isend(buffer22,(lx1+4)*(lx3)*lhalo,mpi_realprec, &
!                      idup,tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(2)=tmpreq
!
!    call mpi_irecv(buffer23,(lx1+4)*(lx3)*lhalo,mpi_realprec,idup,&
!                      tag,MPI_COMM_WORLD,tmpreq,ierr)
!    requests(3)=tmpreq
!
!    call mpi_irecv(buffer24,(lx1+4)*(lx3)*lhalo,mpi_realprec,iddown, &
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

end submodule mpihalo
