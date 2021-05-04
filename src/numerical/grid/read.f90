submodule (grid) grid_read

!use mpimod, only : mpi_realprec

implicit none (type, external)

interface ! readgrid_*.f90
module subroutine get_grid3_raw(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
  e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
character(*), intent(in) :: path
integer, intent(in) :: flagswap
class(curvmesh), intent(inout) :: x
real(wp), intent(out), dimension(:,:,:) :: g1all,g2all,g3all,glatall,nullptsall,rall,thetaall,phiall
real(wp), intent(out), dimension(:,:) :: Incall
real(wp), intent(out), dimension(:,:,:,:) :: e1all,e2all,e3all,erall,ethetaall,ephiall
end subroutine get_grid3_raw

module subroutine get_grid3_hdf5(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
  e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
character(*), intent(in) :: path
integer, intent(in) :: flagswap
class(curvmesh), intent(inout) :: x
real(wp), intent(out), dimension(:,:,:) :: g1all,g2all,g3all,glatall,nullptsall,rall,thetaall,phiall
real(wp), intent(out), dimension(:,:) :: Incall
real(wp), intent(out), dimension(:,:,:,:) :: e1all,e2all,e3all,erall,ethetaall,ephiall
end subroutine get_grid3_hdf5

module subroutine get_grid3_nc4(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
  e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
character(*), intent(in) :: path
integer, intent(in) :: flagswap
class(curvmesh), intent(inout) :: x
real(wp), intent(out), dimension(:,:,:) :: g1all,g2all,g3all,glatall,nullptsall,rall,thetaall,phiall
real(wp), intent(out), dimension(:,:) :: Incall
real(wp), intent(out), dimension(:,:,:,:) :: e1all,e2all,e3all,erall,ethetaall,ephiall
end subroutine get_grid3_nc4
end interface

contains

subroutine get_grid3(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
  e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
!! switchyard file format
character(*), intent(in) :: path
integer, intent(in) :: flagswap
type(curvmesh), intent(inout) :: x
real(wp), intent(out), dimension(:,:,:) :: g1all,g2all,g3all,glatall,nullptsall,rall,thetaall,phiall
real(wp), intent(out), dimension(:,:) :: Incall
real(wp), intent(out), dimension(:,:,:,:) :: e1all,e2all,e3all,erall,ethetaall,ephiall

character(:), allocatable :: fmt

fmt = path(index(path, '.', back=.true.) : len(path))

select case (fmt)
case ('.dat')
  call get_grid3_raw(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
    e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
case ('.h5')
  call get_grid3_hdf5(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
    e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
case ('.nc')
  call get_grid3_nc4(path, flagswap, x, g1all,g2all,g3all,glatall,Incall,nullptsall,&
    e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)
case default
  write(stderr,*) 'grid:read:get_grid3: unknown grid format: ' // fmt
  error stop 2
end select

end subroutine get_grid3


!> Read in the grid information and prep grid object.  Note that there are also module-scope variables
!   that are (redundantly, for convenience) defined based on the grid structure and this procedure 
!   must also set those.  
module procedure read_grid
! subroutine read_grid(indatsize,indatgrid,flagperiodic,x)

!> Read file information for grid; here just the coordinate arrays
!   we also need some information from the input file about what type
!   of grid (e.g. dipole or cartesian) is to be generated.  These routines
!   will also set grid-level array sizes.
  if(mpi_cfg%myid==0) then
    call read_grid_root(indatsize,indatgrid,x1,x2,x3,x2all,x3all)
  else
    call read_grid_workers(x1,x2,x3,x2all,x3all)
  end if
  
  !> Declare grid type that we are dealing with; note lack of matching deallocate assume
  !   that the compiler will deal with it automatically
  select case (gridflag)
    case(0)    ! cartesian
      allocate(cartmesh::x) 
    case(1)    ! dipole
      allocate(dipolemesh::x)
    case default
      error stop ' invalid mesh type specified; cannot instantiate!'
  end select

  !> Create the grid object
  call x%set_coords(x1,x2,x3,x2all,x3all)    ! store coordinate arrays
  call x%init()                              ! allocate space for subgrid variables
  call x%make()                              ! fill auxiliary arrays
  
  !> We need to collect the info for root's fullgrid variables
  if (mpi_cfg%myid==0) then
    call x%init_storage_root()                ! now we have space in type to store full-grid arrays
    call gather_grid_root(x%h1,x%h2,x%h3, &
                          x%h1x1i,x%h2x1i,x%h3x1i, &
                          x%h1x2i,x%h2x2i,x%h3x2i, &
                          x%h1x2i,x%h2x3i,x%h3x3i, &
                          x%r,x%theta,x%phi, &
                          x%alt,x%Bmag,x%glon &
                          x%h1all,x%h2all,x%h3all, &
                          x%h1x1iall,x%h2x1iall,x%h3x1iall, &
                          x%h1x2iall,x%h2x2iall,x%h3x2iall, &
                          x%h1x2iall,x%h2x3iall,x%h3x3iall, &
                          x%rall,x%thetaall,x%phiall, &
                          x%altall,x%Bmagall,x%glonall)     
    !! note that we can fill arrays manually with our own routines rather than use x%set_root, saves temp arrays and memory
  else
    !! gather
    call gather_grid_workers(x%h1,x%h2,x%h3, &
                          x%h1x1i,x%h2x1i,x%h3x1i, &
                          x%h1x2i,x%h2x2i,x%h3x2i, &
                          x%h1x2i,x%h2x3i,x%h3x3i, &
                          x%r,x%theta,x%phi, &
                          x%alt,x%Bmag,x%glon)
  end if
  
  !> Assign periodic or not based on user input
  call x%set_periodic(flagperiodic)
  
  !> Set flags for module scope vars.
  gridflag=x%gridflag
  
  !> Set gravitational fields for module scope vars., use pointers to avoid duplicating data
  g1=>x%g1; g2=>x%g2; g3=>x%g3
  
  !> Make sure we have a sensible x2,3 decomposition of grid
  !> and that parameters aren't impossible
  if(mpi_cfg%myid == 0) call grid_check(x)
end procedure read_grid


!> Read in basic coordinates and grid type (if present) from input file
subroutine read_grid_root(indatsize,indatgrid,x)

!------------------------------------------------------------
!--------SOME CODE DUPLICATION WITH WORKER VERSION - CAN WE
!--------CREATE A COMMON SUBROUTINE THAT ALLOCATES SOME VARS?
!--------THE GRID SIZES MUST ALREADY BE DEFINED IN THE MODULE
!--------, PROBABLY FROM CALLING GRID_SIZE
!------------------------------------------------------------

  character(*), intent(in) :: indatsize,indatgrid
  class(curvmesh), intent(inout) :: x    !does this need to be inout?  I think if anything is unallocated, it does...
  integer u,ierr,iid,ix1,ix2,ix3,icount,icomp!,itell
  
  !print*, 'Entering read_grid_root', lx1,lx2all,lx3all,mpi_cfg%lid2,mpi_cfg%lid3,lx2all/mpi_cfg%lid2
  
  
  !DETERMINE NUMBER OF SLABS AND CORRESPONDING SIZE FOR EACH WORKER
  !NOTE THAT WE WILL ASSUME THAT THE GRID SIZE IS DIVISIBLE BY NUMBER OF WORKERS AS THIS HAS BEEN CHECKED A FEW LINES BELOW
  x%lx1 = lx1
  x%lx2all = lx2all
  x%lx3all = lx3all
  
  
  !> ADJUST THE SIZES OF THE VARIABLES IF LX3ALL==1, SO THAT THE ALLOCATIONS ARE THE CORRECT SIZE
  if (lx3all==1) then
    print *, 'read_grid_root: 2D run: **SWAP** x2 and x3 dims'
    lx3all=lx2all; x%lx3all=lx2all;
    lx2=1
    lx2all=1; x%lx2all=1;
    lx3=lx3all/mpi_cfg%lid
    !! use total number of processes since we only divide in one direction here...
    flagswap=1
  else
    !! non-swapped axes
    if(lx2all==1) then
      print *, 'read_grid_root: 2D run: do not swap x2 and x3 dims.'
      lx2 = 1
      lx3 = lx3all/mpi_cfg%lid
    else
      print *, 'read_grid_root: 3D run'
      !! should divide evenly if generated from process_grid
      lx2 = lx2all/mpi_cfg%lid2
      lx3 = lx3all/mpi_cfg%lid3
    endif
    flagswap=0
  end if
  x%lx2=lx2; x%lx3=lx3
  print '(A,3I6)', 'Grid slab size:  ',lx1,lx2,lx3
  
  if(lx2all > 1 .and. lx3all > 1) then
    if(lx2 == 1 .or. lx3 == 1) error stop "read_grid_root: 3D grids cannot be partitioned with a single MPI image on an axis"
  endif
  flagswap=0
end if
x%lx2=lx2; x%lx3=lx3
print '(A,3I6)', 'Grid slab size:  ',lx1,lx2,lx3

if(lx2all > 1 .and. lx3all > 1) then
  if(lx2 == 1 .or. lx3 == 1) error stop "read_grid_root: 3D grids cannot be partitioned with a single MPI image on an axis"
endif

!> COMMUNICATE THE GRID SIZE TO THE WORKERS SO THAT THEY CAN ALLOCATE SPACE
ierr=0
!! if this is a root-only simulation we don't want to error out
do iid=1,mpi_cfg%lid-1
  call mpi_send(lx2,1,MPI_INTEGER,iid,tag%lx2,MPI_COMM_WORLD,ierr)
  !! need to also pass the lx2all size to all workers to they know
  !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send lx2'

  call mpi_send(lx3,1,MPI_INTEGER,iid,tag%lx3,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send lx3'
end do

if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send grid size'

!TELL WORKERS IF WE'VE SWAPPED DIMENSIONS
ierr=0
do iid=1,mpi_cfg%lid-1
  call mpi_send(flagswap,1,MPI_INTEGER,iid,tag%swap,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send flagswap'
end do

if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send flagswap'


!ALLOCATE SPACE FOR ROOTS SLAB OF DATA
allocate(x%x1(-1:lx1+2))
allocate(x%dx1i(lx1),x%x1i(lx1+1),x%dx1(0:lx1+2))

!FULL GRID X2 VARIABLE
allocate(x%x2all(-1:lx2all+2))
allocate(x%dx2all(0:lx2all+2))
allocate(x%dx2iall(lx2all),x%x2iall(lx2all+1))

!FULL-GRID X3-VARIABLE
allocate(x%x3all(-1:lx3all+2))
allocate(x%dx3all(0:lx3all+2))
allocate(x%x3iall(lx3all+1),x%dx3iall(lx3all))

allocate(x%h1all(-1:lx1+2,-1:lx2all+2,-1:lx3all+2),x%h2all(-1:lx1+2,-1:lx2all+2,-1:lx3all+2), &
         x%h3all(-1:lx1+2,-1:lx2all+2,-1:lx3all+2))    !do we need the ghost cell values?  Yes the divergence in compression term is computed using one ghost cell in each direction!!!!
allocate(x%h1x1iall(1:lx1+1,1:lx2all,1:lx3all),x%h2x1iall(1:lx1+1,1:lx2all,1:lx3all), &
           x%h3x1iall(1:lx1+1,1:lx2all,1:lx3all))
allocate(x%h1x2iall(1:lx1,1:lx2all+1,1:lx3all),x%h2x2iall(1:lx1,1:lx2all+1,1:lx3all), &
           x%h3x2iall(1:lx1,1:lx2all+1,1:lx3all))
allocate(x%h1x3iall(1:lx1,1:lx2all,1:lx3all+1),x%h2x3iall(1:lx1,1:lx2all,1:lx3all+1), &
           x%h3x3iall(1:lx1,1:lx2all,1:lx3all+1))

allocate(x%rall(1:lx1,1:lx2all,1:lx3all),x%thetaall(1:lx1,1:lx2all,1:lx3all),x%phiall(1:lx1,1:lx2all,1:lx3all))
allocate(x%altall(1:lx1,1:lx2all,1:lx3all),x%Bmagall(1:lx1,1:lx2all,1:lx3all))
allocate(x%glonall(1:lx1,1:lx2all,1:lx3all))


!DETERMINE AND ALLOCATE SPACE NEEDED FOR ROOT SUBGRIDS (WORKERS USE SIMILAR ALLOCATE STATEMENTS)
allocate(x%x3(-1:lx3+2))
allocate(x%dx3i(lx3),x%x3i(lx3+1),x%dx3(0:lx3+2))

allocate(x%x2(-1:lx2+2))
allocate(x%dx2i(lx2),x%x2i(lx2+1),x%dx2(0:lx2+2))

allocate(x%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h2(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
!    allocate(x%h1(1:lx1,1:lx2,1:lx3),x%h2(1:lx1,1:lx2,1:lx3),x%h3(1:lx1,1:lx2,1:lx3))
allocate(x%h1x1i(1:lx1+1,1:lx2,1:lx3),x%h2x1i(1:lx1+1,1:lx2,1:lx3),x%h3x1i(1:lx1+1,1:lx2,1:lx3))
allocate(x%h1x2i(1:lx1,1:lx2+1,1:lx3),x%h2x2i(1:lx1,1:lx2+1,1:lx3),x%h3x2i(1:lx1,1:lx2+1,1:lx3))
allocate(x%h1x3i(1:lx1,1:lx2,1:lx3+1),x%h2x3i(1:lx1,1:lx2,1:lx3+1),x%h3x3i(1:lx1,1:lx2,1:lx3+1))

allocate(x%glat(1:lx1,1:lx2,1:lx3),x%glon(1:lx1,1:lx2,1:lx3),x%alt(1:lx1,1:lx2,1:lx3))
allocate(x%r(1:lx1,1:lx2,1:lx3),x%theta(1:lx1,1:lx2,1:lx3),x%phi(1:lx1,1:lx2,1:lx3))

allocate(x%Bmag(1:lx1,1:lx2,1:lx3))
allocate(x%I(1:lx2,1:lx3))
allocate(x%nullpts(1:lx1,1:lx2,1:lx3))

allocate(x%e1(1:lx1,1:lx2,1:lx3,1:3),x%e2(1:lx1,1:lx2,1:lx3,1:3),x%e3(1:lx1,1:lx2,1:lx3,1:3))
allocate(x%er(1:lx1,1:lx2,1:lx3,1:3),x%etheta(1:lx1,1:lx2,1:lx3,1:3),x%ephi(1:lx1,1:lx2,1:lx3,1:3))


!NOW WE NEED TO READ IN THE FULL GRID DATA AND PUT IT INTO THE STRUCTURE.
!IF WE HAVE DONE ANY DIMENSION SWAPPING HERE WE NEED TO TAKE THAT INTO ACCOUNT IN THE VARIABLES THAT ARE BEING READ IN
print *, 'read_grid_root: Starting grid input from file: ',indatgrid

allocate(g1all(lx1,lx2all,lx3all), g2all(lx1,lx2all,lx3all), g3all(lx1,lx2all,lx3all))
!allocate(altall(lx1,lx2all,lx3all), glatall(lx1,lx2all,lx3all), glonall(lx1,lx2all,lx3all))
allocate(glatall(lx1,lx2all,lx3all))
!allocate(Bmagall(lx1,lx2all,lx3all))
allocate(Incall(lx2all,lx3all))
allocate(nullptsall(lx1,lx2all,lx3all))
allocate(e1all(lx1,lx2all,lx3all,3), e2all(lx1,lx2all,lx3all,3), e3all(lx1,lx2all,lx3all,3))
allocate(erall(lx1,lx2all,lx3all,3), ethetaall(lx1,lx2all,lx3all,3), ephiall(lx1,lx2all,lx3all,3))
allocate(rall(lx1,lx2all,lx3all), thetaall(lx1,lx2all,lx3all), phiall(lx1,lx2all,lx3all))

if (flagswap/=1) then
  !! normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
else
  !! this is apparently a 2D grid, so the x2 and x3 dimensions have been/need to be swapped
  !! MZ - may need to change lx2-->lx2all???
  !! Note there that the fortran arrays are the correct size, but the input data are not!!!
  !! This means tmp variable and permutes...

  print *, 'read_grid_root: 2D grid: **PERMUTE** x2 and x3 dimensions'
end if

call get_grid3(indatgrid,flagswap,x,g1all,g2all,g3all,glatall,Incall,nullptsall,&
  e1all,e2all,e3all,erall,ethetaall,ephiall,rall,thetaall,phiall)   !altall and Bmagall now in grid structure (if root)

!! ALLOCATE SPACE FOR ROOTS SUBGRID GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))


!! SEND FULL X1 AND X2 GRIDS TO EACH WORKER (ONLY X3-DIM. IS INVOLVED IN THE MPI
print*, 'read_grid_root: Exchanging grid spacings...'
do iid=1,mpi_cfg%lid-1
  call mpi_send(x%x1,lx1+4,mpi_realprec,iid,tag%x1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x1'
  call mpi_send(x%x2all,lx2all+4,mpi_realprec,iid,tag%x2all,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x2all'
  call mpi_send(x%x3all,lx3all+4,mpi_realprec,iid,tag%x3all,MPI_COMM_WORLD,ierr)
  !! workers may need a copy of this, e.g. for boundary conditions
  !if (ierr/=0) error stop 'failed mpi_send x3all'

  call mpi_send(x%dx1,lx1+3,mpi_realprec,iid,tag%x1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send dx1'
  call mpi_send(x%x1i,lx1+1,mpi_realprec,iid,tag%x1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x1i'
  call mpi_send(x%dx1i,lx1,mpi_realprec,iid,tag%x1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send dx1i'
end do

if (ierr/=0) error stop 'failed mpi_send x grid'

!! NOW SEND THE INFO THAT DEPENDS ON X3 SLAB SIZE
print*, 'read_grid_root: Computing subdomain spacing...'
call bcast_send1D_3(x%x3all,tag%x3,x%x3)
x%dx3 = x%x3(0:lx3+2)-x%x3(-1:lx3+1)
!! computing these avoids extra message passing (could be done for other coordinates, as well)
x%x3i(1:lx3+1)=0.5*(x%x3(0:lx3)+x%x3(1:lx3+1))
x%dx3i=x%x3i(2:lx3+1)-x%x3i(1:lx3)


call bcast_send1D_2(x%x2all,tag%x2,x%x2)
x%dx2 = x%x2(0:lx2+2)-x%x2(-1:lx2+1)
!! computing these avoids extra message passing (could be done for other coordinates
x%x2i(1:lx2+1) = 0.5*(x%x2(0:lx2)+x%x2(1:lx2+1))
x%dx2i=x%x2i(2:lx2+1)-x%x2i(1:lx2)


print*, 'read_grid_root: Dealing with metric factors...'
allocate(mpisendbuf(-1:lx1+2,-1:lx2all+2,-1:lx3all+2),mpirecvbuf(-1:lx1+2,-1:lx2+2,-1:lx3+2))

mpisendbuf=x%h1all
!! since metric factors are pointers they are not guaranteed to be contiguous in memory so pack them into a buffer that is...
call bcast_send3D_ghost(mpisendbuf,tag%h1,mpirecvbuf)
!! special broadcast subroutine to handle 3D arrays with ghost cells
x%h1=mpirecvbuf       !store roots slab of metric factors in its grid structure
mpisendbuf=x%h2all
call bcast_send3D_ghost(mpisendbuf,tag%h2,mpirecvbuf)
x%h2=mpirecvbuf
mpisendbuf=x%h3all
call bcast_send3D_ghost(mpisendbuf,tag%h3,mpirecvbuf)
x%h3=mpirecvbuf

deallocate(mpisendbuf,mpirecvbuf)    !we need different sized buffers below

call bcast_send(x%h1x1iall,tag%h1,x%h1x1i)
!! do the weird sizes here (ie. lx1+1) give issues with MPI?  probably...  No because bcast reads the size off of the variable...
call bcast_send(x%h2x1iall,tag%h2,x%h2x1i)
call bcast_send(x%h3x1iall,tag%h3,x%h3x1i)

call bcast_send3D_x2i(x%h1x2iall,tag%h1,x%h1x2i)
call bcast_send3D_x2i(x%h2x2iall,tag%h2,x%h2x2i)
call bcast_send3D_x2i(x%h3x2iall,tag%h3,x%h3x2i)

call bcast_send3D_x3i(x%h1x3iall,tag%h1,x%h1x3i)
call bcast_send3D_x3i(x%h2x3iall,tag%h2,x%h2x3i)
call bcast_send3D_x3i(x%h3x3iall,tag%h3,x%h3x3i)


print *, 'read_grid_root: Sending gravity, etc...'
call bcast_send(g1all,tag%h1,g1)
call bcast_send(g2all,tag%h2,g2)
call bcast_send(g3all,tag%h3,g3)

call bcast_send(glatall,tag%glat,x%glat)
call bcast_send(x%glonall,tag%glon,x%glon)
call bcast_send(x%altall,tag%alt,x%alt)

call bcast_send(x%Bmagall,tag%Bmag,x%Bmag)
call bcast_send(Incall,tag%inc,x%I)
call bcast_send(nullptsall,tag%null,x%nullpts)

print *, 'read_grid_root: Now sending unit vectors...'
allocate(mpisendbuf(1:lx1,1:lx2all,1:lx3all),mpirecvbuf(1:lx1,1:lx2,1:lx3))
!! why is buffering used/needed here???
do icomp=1,3
  mpisendbuf=e1all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%eunit1,mpirecvbuf)
  x%e1(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e2all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%eunit2,mpirecvbuf)
  x%e2(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e3all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%eunit3,mpirecvbuf)
  x%e3(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=erall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%er,mpirecvbuf)
  x%er(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ethetaall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%etheta,mpirecvbuf)
  x%etheta(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ephiall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tag%ephi,mpirecvbuf)
  x%ephi(:,:,:,icomp)=mpirecvbuf
end do
deallocate(mpisendbuf,mpirecvbuf)

call bcast_send(rall,tag%r,x%r)
call bcast_send(thetaall,tag%theta,x%theta)
call bcast_send(phiall,tag%phi,x%phi)
print *,  'read_grid_root: Done sending slabbed variables to workers...'


!> COUNT THE NUMBER OF NULL GRID POINTS AND GENERATE A LIST OF NULL INDICES FOR LATER USE
x%lnull=0;
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3) > 0.5_wp) then
        x%lnull=x%lnull+1
      end if
    end do
  end do
  if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send grid size'
  
  !TELL WORKERS IF WE'VE SWAPPED DIMENSIONS
  ierr=0
  do iid=1,mpi_cfg%lid-1
    call mpi_send(flagswap,1,MPI_INTEGER,iid,tag%swap,MPI_COMM_WORLD,ierr)
    !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send flagswap'
  end do
end do
print *, 'read_grid_root: Done computing null grid points...  Process:  ',mpi_cfg%myid,' has:  ',x%lnull


!COMPUTE DIFFERENTIAL DISTANCES ALONG EACH DIRECTION (TO BE USED IN TIME STEP DETERMINATION...
block
real(wp), dimension(1:lx1,1:lx2,1:lx3) :: tmpdx
allocate(x%dl1i(1:lx1,1:lx2,1:lx3),x%dl2i(1:lx1,1:lx2,1:lx3),x%dl3i(1:lx1,1:lx2,1:lx3))

tmpdx=spread(spread(x%dx1i,2,lx2),3,lx3)
x%dl1i=tmpdx*x%h1(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx2i,1,lx1),3,lx3)
x%dl2i=tmpdx*x%h2(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx3i,1,lx1),2,lx2)
x%dl3i=tmpdx*x%h3(1:lx1,1:lx2,1:lx3)
end block

!    !THIS CODE BLOCK HAS ROOT "PARROT" THE GRID TO A FILE FOR DEBUGGING PURPOSES.
!    !THIS BLOCK MUST BE HERE DUE TO THE FACT MANY OF THE ALL-GRID VARIABLES ARE ONLY
!    !IN SCOPE INSIDE THIS FUNCTION BEFORE THE DEALLOCATE STATEMENT BELOW
!    open(newunit=outunit,file='testsize.dat',status='replace',form='unformatted', &
!         access='stream', action='write')
!    write(outunit) lx1,lx2,lx3all    !note that these sizes exxclude ghost cells
!    close(outunit)
!
!    open(outunit,file='testgrid.dat',status='replace',form='unformatted', &
!         access='stream')
!    write(outunit) x%x1,x%x1i,x%dx1,x%dx1i    !I guess the number of elements is obtained from the size of the arrays???
!    write(outunit) x%x2,x%x2i,x%dx2,x%dx2i
!    write(outunit) x%x3all,x%x3iall,x%dx3all,x%dx3iall
!    write(outunit) x%h1all,x%h2all,x%h3all
!    write(outunit) x%h1x1iall,x%h2x1iall,x%h3x1iall
!    write(outunit) x%h1x2iall,x%h2x2iall,x%h3x2iall
!    write(outunit) x%h1x3iall,x%h2x3iall,x%h3x3iall
!    write(outunit) g1all,g2all,g3all
!    write(outunit) altall,glatall,glonall
!    write(outunit) Bmagall
!    write(outunit) Incall
!    write(outunit) nullptsall
!    close(outunit)
!    print *, 'Done creating copy of grid...'
!


!    !FLAG THE GRID AS PERIODIC, IF REQUESTED; IF SO PERIODICITY WILL BE ASSUMED IN THE X3-DIRECTION
!    if (flagperiodic==1) then
!      x%flagper=1
!    else
!      x%flagper=0
!    end if
!

!> DEALLOCATE ANY FULL GRID VARIABLES THAT ARE NO LONGER NEEDED
deallocate(g1all,g2all,g3all,glatall,Incall,nullptsall,rall,thetaall,phiall)
!! altall and Bmagall now in x structure

end subroutine read_grid_root


subroutine read_grid_workers(x)

class(curvmesh), intent(inout) :: x
integer :: ix1,ix2,ix3,icount,icomp, ierr

!GET ROOTS MESSAGE WITH THE SIZE OF THE GRID WE ARE TO RECEIVE
call mpi_recv(lx2,1,MPI_INTEGER,0,tag%lx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv lx2'
call mpi_recv(lx3,1,MPI_INTEGER,0,tag%lx3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv lx3'

x%lx1=lx1; x%lx2=lx2; x%lx2all=lx2all; x%lx3all=lx3all; x%lx3=lx3

!ROOT NEEDS TO TELL US WHETHER WE'VE SWAPPED DIMENSIONS SINCE THIS AFFECTS HOW CURRENTS ARE COMPUTED
call mpi_recv(flagswap,1,MPI_INTEGER,0,tag%swap,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv flagswap'

if (flagswap==1) then
  x%lx3all=lx2all
  x%lx2all=lx3all
  lx2all=x%lx2all
  lx3all=x%lx3all
end if

!ALLOCATE SPACE FOR WORKER'S GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))

end subroutine read_grid_workers


!> pull full grid vars. from workers into root arrays
subroutine gather_grid_root(h1,h2,h3, &
                        h1x1i,h2x1i,h3x1i, &
                        h1x2i,h2x2i,h3x2i, &
                        h1x2i,h2x3i,h3x3i, &
                        r,theta,phi, &
                        alt,Bmag,glon &
                        h1all,h2all,h3all, &
                        h1x1iall,h2x1iall,h3x1iall, &
                        h1x2iall,h2x2iall,h3x2iall, &
                        h1x2iall,h2x3iall,h3x3iall, &
                        rall,thetaall,phiall, &
                        altall,Bmagall,glonall)
  real(wp), dimension(:,:,:), intent(in) :: h1,h2,h3
  real(wp), dimension(:,:,:), intent(in) :: h1x1i,h2x1i,h3x1i
  real(wp), dimension(:,:,:), intent(in) :: h1x2i,h2x2i,h3x2i
  real(wp), dimension(:,:,:), intent(in) :: h1x3i,h2x3i,h3x3i
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(:,:,:), intent(in) :: alt,Bmag,glon
  real(wp), dimension(:,:,:), intent(out) :: h1all,h2all,h3all
  real(wp), dimension(:,:,:), intent(out) :: h1x1iall,h2x1iall,h3x1iall
  real(wp), dimension(:,:,:), intent(out) :: h1x2iall,h2x2iall,h3x2iall
  real(wp), dimension(:,:,:), intent(out) :: h1x3iall,h2x3iall,h3x3iall
  real(wp), dimension(:,:,:), intent(out) :: rall,thetaall,phiall
  real(wp), dimension(:,:,:), intent(out) :: altall,Bmagall,glonall

  call gather_recv3D_ghost(h1,tag%h1,h1all)
  call gather_recv3D_ghost(h2,tag%h2,h2all)
  call gather_recv3D_ghost(h3,tag%h3,h3all)

  call gather_recv(h1x1i,tag%h1,h1x1iall)
  call gather_recv(h2x1i,tag%h2,h2x1iall)
  call gather_recv(h3x1i,tag%h3,h3x1iall)

  call gather_recv3D_x2i(h1x2i,tag%h1,h1x2iall)
  call gather_recv3D_x2i(h2x2i,tag%h2,h2x2iall)
  call gather_recv3D_x2i(h3x2i,tag%h3,h3x2iall)

  call gather_recv3D_x3i(h1x3i,tag%h1,h1x3iall)
  call gather_recv3D_x3i(h2x3i,tag%h2,h2x3iall)
  call gather_recv3D_x3i(h3x3i,tag%h3,h3x3iall)

  call gather_recv(r,tag%r,rall)
  call gather_recv(theta,tag%theta,thetaall)
  call gather_recv(phi,tag%phi,phiall)

  call gather_recv(alt,tag%alt,altall)
  call gather_recv(Bmag,tag%Bmag,Bmagall)
  call gather_recv(glon,tag%glon,glonall)
end subroutien gather_grid_root


!> send full grid vars. to root
subroutine gather_grid_workers(h1,h2,h3, &
                        h1x1i,h2x1i,h3x1i, &
                        h1x2i,h2x2i,h3x2i, &
                        h1x2i,h2x3i,h3x3i, &
                        r,theta,phi, &
                        alt,Bmag,glon)
  real(wp), dimension(:,:,:), intent(in) :: h1,h2,h3
  real(wp), dimension(:,:,:), intent(in) :: h1x1i,h2x1i,h3x1i
  real(wp), dimension(:,:,:), intent(in) :: h1x2i,h2x2i,h3x2i
  real(wp), dimension(:,:,:), intent(in) :: h1x3i,h2x3i,h3x3i
  real(wp), dimension(:,:,:), intent(in) :: r,theta,phi
  real(wp), dimension(:,:,:), intent(in) :: alt,Bmag,glon

  call gather_send3D_ghost(h1,tag%h1)
  call gather_send3D_ghost(h2,tag%h2)
  call gather_send3D_ghost(h3,tag%h3)

  call gather_send(h1x1i,tag%h1)
  call gather_send(h2x1i,tag%h2)
  call gather_send(h3x1i,tag%h3)

  call gather_send3D_x2i(h1x2i,tag%h1)
  call gather_send3D_x2i(h2x2i,tag%h2)
  call gather_send3D_x2i(h3x2i,tag%h3)

  call gather_send3D_x3i(h1x3i,tag%h1)
  call gather_send3D_x3i(h2x3i,tag%h2)
  call gather_send3D_x3i(h3x3i,tag%h3)

  call gather_send(r,tag%r)
  call gather_send(theta,tag%theta)
  call gather_send(phi,tag%phi)

  call gather_send(alt,tag%alt)
  call gather_send(Bmag,tag%Bmag)
  call gather_send(glon,tag%glon)
end subroutine gather_grid_workers

end submodule grid_read
