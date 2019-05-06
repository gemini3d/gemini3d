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
