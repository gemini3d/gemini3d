submodule (grid) grid_read
implicit none
contains


module procedure read_grid
! subroutine read_grid(indatsize,indatgrid,flagperiodic,x)
!! READS GRID INFORMATION FROM A BBINARY FILE AND ALLOCATES
!! GRID STRUCTURE.  NOTE THAT THE INPUT FILE DOES NOT,
!! BY DEFAULT INCLUDE THE GHOST CELLS.  THIS FUNCTION
!! IS A WRAPPER WHICH PASSES WORK OFF TO APPROPRIATE
!! SUBROUTINE DEPENDING ON WHETHER IT IS CALLED BY ROOT
!! OR NOT.   THIS SUBRoutINE ALSO CLASSIFIES THE GRID.



!> READ IN THE GRID DATA
if(myid==0) then
  call read_grid_root(indatsize,indatgrid,x)
else
  call read_grid_workers(x)
end if


!FLAG THE GRID AS PERIODIC, IF REQUESTED; IF SO PERIODICITY WILL BE ASSUMED
!IN THE X3-DIRECTION, NOTE BOTH ROOT AND WORKERS MUST DO THIS!!!
if (flagperiodic==1) then
  x%flagper=.true.
else
  x%flagper=.false.
end if


!DETERMINE THE TYPE OF GRID WE HAVE AND SET AN APPROPRIATE FLAG
!FIXME:  this needs to be done once, globally and also needs to be more
!robust...
if (abs(x%alt(1,1,1)-x%alt(lx1,1,1))<100d3) then    !closed dipole grid
  gridflag=0
else if (x%alt(1,1,1)>x%alt(2,1,1)) then    !open dipole grid with inverted structure wrt altitude
  gridflag=1
else    !something different (viz. non-inverted - lowest altitudes at the logical bottom of the grid)
  gridflag=2
end if

end procedure read_grid


subroutine read_grid_root(indatsize,indatgrid,x)

!------------------------------------------------------------
!--------SOME CODE DUPLICATION WITH WORKER VERSION - CAN WE
!--------CREATE A COMMON SUBROUTINE THAT ALLOCATES SOME VARS?
!--------THE GRID SIZES MUST ALREADY BE DEFINED IN TEH MODULE
!--------, PROBABLY FROM CALLING GRID_SIZE
!------------------------------------------------------------

character(*), intent(in) :: indatsize,indatgrid
type(curvmesh), intent(inout) :: x    !does this need to be inout?  I think if anything is unallocated, it does...

integer u,ierr,iid,ix1,ix2,ix3,icount,icomp!,itell

!NOTE THAT HAVING THESE ARE LOCAL (TEMPORARY) VARS. PREVENTS ROOT FROM WRITING ENTIRE GRID TO FILE AT SOME LATER POINT...
real(wp), dimension(:,:,:), allocatable :: g1all,g2all,g3all   !to temporarily store input data to be distributed
real(wp), dimension(:,:,:), allocatable :: altall,glatall,glonall
real(wp), dimension(:,:,:), allocatable :: rall,thetaall,phiall
real(wp), dimension(:,:,:), allocatable :: Bmagall
real(wp), dimension(:,:), allocatable :: Incall
real(wp), dimension(:,:,:), allocatable :: nullptsall
real(wp), dimension(:,:,:,:), allocatable :: e1all,e2all,e3all,erall,ethetaall,ephiall    !might be best to have a tmp vector array...
real(wp), dimension(:,:,:), allocatable :: mpisendbuf
real(wp), dimension(:,:,:), allocatable :: mpirecvbuf

!print*, 'Entering read_grid_root', lx1,lx2all,lx3all,lid2,lid3,lx2all/lid2


!DETERMINE NUMBER OF SLABS AND CORRESPONDING SIZE FOR EACH WORKER
!NOTE THAT WE WILL ASSUME THAT THE GRID SIZE IS DIVISIBLE BY NUMBER OF WORKERS AS THIS HAS BEEN CHECKED A FEW LINES BELOW
x%lx1=lx1; x%lx2all=lx2all; x%lx3all=lx3all;


!> ADJUST THE SIZES OF THE VARIABLES IF LX3ALL==1, SO THAT THE ALLOCATIONS ARE THE CORRECT SIZE
if (lx3all==1) then
  print *, 'Detected a 2D run request...  swapping x2 and x3 sizes to maintain parallelization.'
  lx3all=lx2all; x%lx3all=lx2all;
  lx2=1
  lx2all=1; x%lx2all=1;
  lx3=lx3all/lid
  !! use total number of processes since we only divide in one direction here...
  flagswap=1
else
  if(lx2all==1) then
    print *, 'Detected a 2D run request...  ***NOT*** swapping x2 and x3 sizes to maintain parallelization.'
    lx2=1
    lx3=lx3all/lid
  else
    print *, 'Detected a 3D run requuest...'
    lx2=lx2all/lid2
    !! should divide evenly if generated from mpigrid
    lx3=lx3all/lid3
  endif
  flagswap=0
end if
x%lx2=lx2; x%lx3=lx3
print *, 'Grid slab size:  ',lx1,lx2,lx3


!> COMMUNICATE THE GRID SIZE TO THE WORKERS SO THAT THEY CAN ALLOCATE SPACE
ierr=0
!! if this is a root-only simulation we don't want to error out
do iid=1,lid-1
  call mpi_send(lx2,1,MPI_INTEGER,iid,taglx2,MPI_COMM_WORLD,ierr)
  !! need to also pass the lx2all size to all workers to they know
  !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send lx2'

  call mpi_send(lx3,1,MPI_INTEGER,iid,taglx3,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send lx3'
end do

if (ierr/=0) error stop 'grid:read_grid_root failed mpi_send grid size'

!TELL WORKERS IF WE'VE SWAPPED DIMENSIONS
ierr=0
do iid=1,lid-1
  call mpi_send(flagswap,1,MPI_INTEGER,iid,tagswap,MPI_COMM_WORLD,ierr)
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
print *, 'Starting grid input from file: ',indatgrid
open(newunit=u,file=indatgrid,status='old',form='unformatted',access='stream', action='read')
if (flagswap/=1) then                     !normal (i.e. full 3D) grid ordering, or a 2D grid with 1 element naturally in the second dimension
  read(u) x%x1,x%x1i,x%dx1,x%dx1i
  read(u) x%x2all,x%x2iall,x%dx2all,x%dx2iall
  read(u) x%x3all,x%x3iall,x%dx3all,x%dx3iall
  read(u) x%h1all,x%h2all,x%h3all
  read(u) x%h1x1iall,x%h2x1iall,x%h3x1iall
  read(u) x%h1x2iall,x%h2x2iall,x%h3x2iall
  read(u) x%h1x3iall,x%h2x3iall,x%h3x3iall

  allocate(g1all(lx1,lx2all,lx3all),g2all(lx1,lx2all,lx3all),g3all(lx1,lx2all,lx3all))
  read(u) g1all,g2all,g3all

  allocate(altall(lx1,lx2all,lx3all),glatall(lx1,lx2all,lx3all),glonall(lx1,lx2all,lx3all))
  read(u) altall,glatall,glonall

  allocate(Bmagall(lx1,lx2all,lx3all))
  read(u) Bmagall

  allocate(Incall(lx2all,lx3all))
  read(u) Incall

  allocate(nullptsall(lx1,lx2all,lx3all))
  read(u) nullptsall

  allocate(e1all(lx1,lx2all,lx3all,3))
  read(u) e1all
  allocate(e2all(lx1,lx2all,lx3all,3))
  read(u) e2all
  allocate(e3all(lx1,lx2all,lx3all,3))
  read(u) e3all

  allocate(erall(lx1,lx2all,lx3all,3))
  read(u) erall
  allocate(ethetaall(lx1,lx2all,lx3all,3))
  read(u) ethetaall
  allocate(ephiall(lx1,lx2all,lx3all,3))
  read(u) ephiall

  allocate(rall(lx1,lx2all,lx3all))
  read(u) rall
  allocate(thetaall(lx1,lx2all,lx3all))
  read(u) thetaall
  allocate(phiall(lx1,lx2all,lx3all))
  read(u) phiall

  x%rall=rall; x%thetaall=thetaall; x%phiall=phiall;
else
  !! this is apparently a 2D grid, so the x2 and x3 dimensions have been/need to be swapped
  !! MZ - may need to change lx2-->lx2all???
  !! Note there that the fortran arrays are the correct size, but the input data are not!!!
  !! This means tmp variable and permutes...

  print *, 'Detected a 2D grid, so will permute the dimensions of the input'

  read(u) x%x1,x%x1i,x%dx1,x%dx1i                !< x1 untouched
  read(u) x%x3all,x%x3iall,x%dx3all,x%dx3iall    !< for a 3D grid this is x2, but now considered x3(all)
  read(u) x%x2all,x%x2iall,x%dx2all,x%dx2iall    !< formerly x3, now x2

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(htmp(-1:lx1+2,-1:lx3all+2,-1:lx2all+2))    !< this stores the input metric factors which are swapped x2/x3 vs. what this simulation will use
  read(u) htmp
  x%h1all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])
  read(u) htmp        !< this would be h3, but with the input structure shape
  x%h3all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])   !< permute the dimensions of the array 3 --> 2, 2 --> 3
  read(u) htmp        !< this would be h3, but with the input structure shape
  x%h2all=reshape(htmp,[lx1+4,lx2all+4,lx3all+4],order=[1,3,2])
  end block

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(htmp(1:lx1+1,1:lx3all,1:lx2all))    !input 2 vs. 3 dimensions swapped from this program
  read(u) htmp
  x%h1x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  x%h3x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  x%h2x1iall=reshape(htmp,[lx1+1,lx2all,lx3all],order=[1,3,2])
  end block

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(htmp(1:lx1,1:lx3all+1,1:lx2all))
  read(u) htmp
  x%h1x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])    !Note also that the x2 interface from teh input file is x3i in this simulation
  read(u) htmp
  x%h3x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
  read(u) htmp
  x%h2x3iall=reshape(htmp,[lx1,lx2all,lx3all+1],order=[1,3,2])
  end block

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(htmp(1:lx1,1:lx3all,1:lx2all+1))
  read(u) htmp
  x%h1x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  read(u) htmp
  x%h3x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  read(u) htmp
  x%h2x2iall=reshape(htmp,[lx1,lx2all+1,lx3all],order=[1,3,2])
  end block

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(g1all(lx1,lx2all,lx3all),g2all(lx1,lx2all,lx3all),g3all(lx1,lx2all,lx3all))
  allocate(htmp(lx1,lx3all,lx2all))
  read(u) htmp
  g1all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  g3all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  g2all=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])

  allocate(altall(lx1,lx2all,lx3all),glatall(lx1,lx2all,lx3all),glonall(lx1,lx2all,lx3all))
  read(u) htmp
  altall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  glatall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  glonall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])

  allocate(Bmagall(lx1,lx2all,lx3all))
  read(u) htmp
  Bmagall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  deallocate(htmp)
  end block

  block
  real(wp), dimension(:,:), allocatable :: htmp2D
  allocate(Incall(lx2all,lx3all))
  allocate(htmp2D(lx3all,lx2all))
  read(u) htmp2D
  Incall=reshape(htmp2D,[lx2all,lx3all],order=[2,1])
  end block


  !inquire(u, pos=itell)
  !print *,'file pos before read',itell
  !read(u) nullptsall
  !print *,'lx1',lx1,'lx2',lx2,'lx3all',lx3all
  !print *,'shape(nullptsall)',shape(nullptsall)
  !inquire(u, pos=itell)
  !print *,'file pos after read',itell
! FIXME BROKEN!
  !allocate(nullptsall(lx1,lx2,lx3all))
  !print *,shape(nullptsall)
  !stop

  ! FIXME would be like this, but this doesn't work.
  !allocate(nullptsall(lx1,lx2,lx3all))
  !read(u) nullptsall

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(htmp(lx1,lx3all,lx2all))
  read(u) htmp
  nullptsall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  deallocate(htmp)
  end block

  block
  real(wp), dimension(:,:,:,:), allocatable :: htmp4D
  allocate(e1all(lx1,lx2all,lx3all,3))
  allocate(htmp4D(lx1,lx3all,lx2all,3))
  read(u) htmp4D
  e1all=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  allocate(e3all(lx1,lx2all,lx3all,3))          !swap the x2/x3 unit vectors
  read(u) htmp4D
  e3all=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  allocate(e2all(lx1,lx2all,lx3all,3))
  read(u) htmp4D
  e2all=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  allocate(erall(lx1,lx2all,lx3all,3))
  read(u) htmp4D
  erall=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  allocate(ethetaall(lx1,lx2all,lx3all,3))
  read(u) htmp4D
  ethetaall=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])

  allocate(ephiall(lx1,lx2all,lx3all,3))
  read(u) htmp4D
  ephiall=reshape(htmp4D,[lx1,lx2all,lx3all,3],order=[1,3,2,4])
  end block

  block
  real(wp), dimension(:,:,:), allocatable :: htmp
  allocate(rall(lx1,lx2all,lx3all),thetaall(lx1,lx2all,lx3all),phiall(lx1,lx2all,lx3all))
  allocate(htmp(lx1,lx3all,lx2all))
  read(u) htmp
  rall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  thetaall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  read(u) htmp
  phiall=reshape(htmp,[lx1,lx2all,lx3all],order=[1,3,2])
  end block

  x%rall=rall; x%thetaall=thetaall; x%phiall=phiall;
end if
close(u)
print *, 'Done reading in grid data...'


!ALLOCATE SPACE FOR ROOTS SUBGRID GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))


!SEND FULL X1 AND X2 GRIDS TO EACH WORKER (ONLY X3-DIM. IS INVOLVED IN THE MPI
print*, 'Exchanging grid spacings...'
do iid=1,lid-1
  call mpi_send(x%x1,lx1+4,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x1'
  call mpi_send(x%x2all,lx2all+4,mpi_realprec,iid,tagx2all,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x2all'
  call mpi_send(x%x3all,lx3all+4,mpi_realprec,iid,tagx3all,MPI_COMM_WORLD,ierr)
  !! workers may need a copy of this, e.g. for boudnary conditions
  !if (ierr/=0) error stop 'failed mpi_send x3all'

  call mpi_send(x%dx1,lx1+3,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send dx1'
  call mpi_send(x%x1i,lx1+1,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send x1i'
  call mpi_send(x%dx1i,lx1,mpi_realprec,iid,tagx1,MPI_COMM_WORLD,ierr)
  !if (ierr/=0) error stop 'failed mpi_send dx1i'
end do

if (ierr/=0) error stop 'failed mpi_send x grid'

!NOW SEND THE INFO THAT DEPENDS ON X3 SLAB SIZE
print*, 'Computing subdomain spacing...'
call bcast_send1D_3(x%x3all,tagx3,x%x3)
x%dx3=x%x3(0:lx3+2)-x%x3(-1:lx3+1)     !computing these avoids extra message passing (could be done for other coordinates, as well)
x%x3i(1:lx3+1)=0.5*(x%x3(0:lx3)+x%x3(1:lx3+1))
x%dx3i=x%x3i(2:lx3+1)-x%x3i(1:lx3)


call bcast_send1D_2(x%x2all,tagx2,x%x2)
x%dx2=x%x2(0:lx2+2)-x%x2(-1:lx2+1)     !computing these avoids extra message passing (could be done for other coordinates
x%x2i(1:lx2+1)=0.5*(x%x2(0:lx2)+x%x2(1:lx2+1))
x%dx2i=x%x2i(2:lx2+1)-x%x2i(1:lx2)


print*, 'Dealing with metric factors...'
allocate(mpisendbuf(-1:lx1+2,-1:lx2all+2,-1:lx3all+2),mpirecvbuf(-1:lx1+2,-1:lx2+2,-1:lx3+2))

mpisendbuf=x%h1all    !since metric factors are pointers they are not gauranteed to be contiguous in memory so pack them into a buffer that is...
call bcast_send3D_ghost(mpisendbuf,tagh1,mpirecvbuf)    !special broadcast subroutine to handle 3D arrays with ghost cells
x%h1=mpirecvbuf       !store roots slab of metric factors in its grid structure
mpisendbuf=x%h2all
call bcast_send3D_ghost(mpisendbuf,tagh2,mpirecvbuf)
x%h2=mpirecvbuf
mpisendbuf=x%h3all
call bcast_send3D_ghost(mpisendbuf,tagh3,mpirecvbuf)
x%h3=mpirecvbuf

deallocate(mpisendbuf,mpirecvbuf)    !we need different sized buffers below

call bcast_send(x%h1x1iall,tagh1,x%h1x1i)    !do the weird sizes here (ie. lx1+1) give issues with MPI?  probably...  No because bcast reads the size off of the variable...
call bcast_send(x%h2x1iall,tagh2,x%h2x1i)
call bcast_send(x%h3x1iall,tagh3,x%h3x1i)

call bcast_send3D_x2i(x%h1x2iall,tagh1,x%h1x2i)
call bcast_send3D_x2i(x%h2x2iall,tagh2,x%h2x2i)
call bcast_send3D_x2i(x%h3x2iall,tagh3,x%h3x2i)

call bcast_send3D_x3i(x%h1x3iall,tagh1,x%h1x3i)
call bcast_send3D_x3i(x%h2x3iall,tagh2,x%h2x3i)
call bcast_send3D_x3i(x%h3x3iall,tagh3,x%h3x3i)


print*, 'Sending gravity, etc...'
call bcast_send(g1all,tagh1,g1)
call bcast_send(g2all,tagh2,g2)
call bcast_send(g3all,tagh3,g3)

call bcast_send(glatall,tagglat,x%glat)
call bcast_send(glonall,tagglon,x%glon)
call bcast_send(altall,tagalt,x%alt)

call bcast_send(Bmagall,tagBmag,x%Bmag)
call bcast_send(Incall,taginc,x%I)
call bcast_send(nullptsall,tagnull,x%nullpts)

print*, 'Now sending unit vectors...'
allocate(mpisendbuf(1:lx1,1:lx2all,1:lx3all),mpirecvbuf(1:lx1,1:lx2,1:lx3))    !why is buffering used/needed here???
do icomp=1,3
  mpisendbuf=e1all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tageunit1,mpirecvbuf)
  x%e1(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e2all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tageunit2,mpirecvbuf)
  x%e2(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=e3all(:,:,:,icomp)
  call bcast_send(mpisendbuf,tageunit3,mpirecvbuf)
  x%e3(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=erall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tager,mpirecvbuf)
  x%er(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ethetaall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tagetheta,mpirecvbuf)
  x%etheta(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  mpisendbuf=ephiall(:,:,:,icomp)
  call bcast_send(mpisendbuf,tagephi,mpirecvbuf)
  x%ephi(:,:,:,icomp)=mpirecvbuf
end do
deallocate(mpisendbuf,mpirecvbuf)

call bcast_send(rall,tagr,x%r)
call bcast_send(thetaall,tagtheta,x%theta)
call bcast_send(phiall,tagphi,x%phi)
print *,  'Done sending slabbed variables to workers...'


!COUNT THE NUMBER OF NULL GRID POINTS AND GENERATE A LIST OF NULL INDICES FOR LATER USE
x%lnull=0;
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%lnull=x%lnull+1
      end if
    end do
  end do
end do
allocate(x%inull(1:x%lnull,1:3))

icount=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%inull(icount,:)=[ix1,ix2,ix3]
        icount=icount+1
      end if
    end do
  end do
end do
print *, 'Done computing null grid points...  Process:  ',myid,' has:  ',x%lnull


!COMPUTE DIFFERENTIAL DISTANCES ALONG EACH DIRECTION (TO BE USED IN TIME STEP DETERMINATION...
block
real(wp), dimension(:,:,:), allocatable :: tmpdx

allocate(x%dl1i(1:lx1,1:lx2,1:lx3),x%dl2i(1:lx1,1:lx2,1:lx3),x%dl3i(1:lx1,1:lx2,1:lx3),tmpdx(1:lx1,1:lx2,1:lx3))

tmpdx=spread(spread(x%dx1i,2,lx2),3,lx3)
x%dl1i=tmpdx*x%h1(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx2i,1,lx1),3,lx3)
x%dl2i=tmpdx*x%h2(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx3i,1,lx1),2,lx2)
x%dl3i=tmpdx*x%h3(1:lx1,1:lx2,1:lx3)
end block

!    !THIS CODE BLOCK HAS ROOT "PARROT" THE GRID TO A FILE FOR DEBUGGING PURPOSES.
!    !THIS BLOCK MUST BE HERE DUE TO THE FACT MANY OF HTE ALL-GRID VARIABLES ARE ONLY
!    !IN SCOPE INSIDE THIS FUNCITON BEFORE THE DEALLOCATE STATEMENT BELOW
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

!DEALLOCATE ANY FULL GRID VARIABLES THAT ARE NO LONGER NEEDED
deallocate(g1all,g2all,g3all,altall,glatall,glonall,Bmagall,Incall,nullptsall,rall,thetaall,phiall)

end subroutine read_grid_root


subroutine read_grid_workers(x)

type(curvmesh), intent(inout) :: x

integer :: ix1,ix2,ix3,icount,icomp, ierr

!GET ROOTS MESSAGE WITH THE SIZE OF THE GRID WE ARE TO RECEIVE
call mpi_recv(lx2,1,MPI_INTEGER,0,taglx2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv lx2'
call mpi_recv(lx3,1,MPI_INTEGER,0,taglx3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv lx3'

x%lx1=lx1; x%lx2=lx2; x%lx2all=lx2all; x%lx3all=lx3all; x%lx3=lx3

!ROOT NEEDS TO TELL US WHETHER WE'VE SWAPPED DIMENSIONS SINCE THIS AFFECTS HOW CURRENTS ARE COMPUTED
call mpi_recv(flagswap,1,MPI_INTEGER,0,tagswap,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv flagswap'

if (flagswap==1) then
  x%lx3all=lx2all
  x%lx2all=lx3all
  lx2all=x%lx2all
  lx3all=x%lx3all
end if


!ALLOCATE SPACE FOR MY SLAB OF DATA
allocate(x%x1(-1:lx1+2))
allocate(x%dx1i(lx1),x%x1i(lx1+1),x%dx1(0:lx1+2))

allocate(x%x2(-1:lx2+2))
allocate(x%dx2i(lx2),x%x2i(lx2+1),x%dx2(0:lx2+2))
allocate(x%x2all(-1:lx2all+2))

!DETERMINE AND ALLOCATE SPACE NEEDED FOR WORKERS SUBGRIDS
allocate(x%x3(-1:lx3+2))
allocate(x%dx3i(lx3),x%x3i(lx3+1),x%dx3(0:lx3+2))
allocate(x%x3all(-1:lx3all+2))

allocate(x%h1(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h2(-1:lx1+2,-1:lx2+2,-1:lx3+2),x%h3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
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


!ALLOCATE SPACE FOR WORKER'S GRAVITATIONAL FIELD
allocate(g1(1:lx1,1:lx2,1:lx3),g2(1:lx1,1:lx2,1:lx3),g3(1:lx1,1:lx2,1:lx3))


!RECEIVE GRID DATA FROM ROOT
call mpi_recv(x%x1,lx1+4,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv x1'
call mpi_recv(x%x2all,lx2all+4,mpi_realprec,0,tagx2all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv x2all'
call mpi_recv(x%x3all,lx3all+4,mpi_realprec,0,tagx3all,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv x3all'
call mpi_recv(x%dx1,lx1+3,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv dx1'
call mpi_recv(x%x1i,lx1+1,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv x1i'
call mpi_recv(x%dx1i,lx1,mpi_realprec,0,tagx1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr/=0) error stop 'failed mpi_recv dx1i'


call bcast_recv1D_3(x%x3,tagx3)
x%dx3=x%x3(0:lx3+2)-x%x3(-1:lx3+1)     !computing these avoids extra message passing (could be done for other coordinates
x%x3i(1:lx3+1)=0.5*(x%x3(0:lx3)+x%x3(1:lx3+1))
x%dx3i=x%x3i(2:lx3+1)-x%x3i(1:lx3)


call bcast_recv1D_2(x%x2,tagx2)
x%dx2=x%x2(0:lx2+2)-x%x2(-1:lx2+1)
x%x2i(1:lx2+1)=0.5*(x%x2(0:lx2)+x%x2(1:lx2+1))
x%dx2i=x%x2i(2:lx2+1)-x%x2i(1:lx2)


block
real(wp), dimension(:,:,:), allocatable :: mpirecvbuf
allocate(mpirecvbuf(-1:lx1+2,-1:lx2+2,-1:lx3+2))

call bcast_recv3D_ghost(mpirecvbuf,tagh1)
x%h1=mpirecvbuf
call bcast_recv3D_ghost(mpirecvbuf,tagh2)
x%h2=mpirecvbuf
call bcast_recv3D_ghost(mpirecvbuf,tagh3)
x%h3=mpirecvbuf
end block

call bcast_recv(x%h1x1i,tagh1)
call bcast_recv(x%h2x1i,tagh2)
call bcast_recv(x%h3x1i,tagh3)

call bcast_recv3D_x2i(x%h1x2i,tagh1)
call bcast_recv3D_x2i(x%h2x2i,tagh2)
call bcast_recv3D_x2i(x%h3x2i,tagh3)

call bcast_recv3D_x3i(x%h1x3i,tagh1)
call bcast_recv3D_x3i(x%h2x3i,tagh2)
call bcast_recv3D_x3i(x%h3x3i,tagh3)

call bcast_recv(g1,tagh1)
call bcast_recv(g2,tagh2)
call bcast_recv(g3,tagh3)

call bcast_recv(x%glat,tagglat)
call bcast_recv(x%glon,tagglon)
call bcast_recv(x%alt,tagalt)

call bcast_recv(x%Bmag,tagBmag)
call bcast_recv(x%I,taginc)           !only time that we need to exchange 2D array data, I think
call bcast_recv(x%nullpts,tagnull)    !note that this is not an integer array

block
real(wp), dimension(:,:,:), allocatable :: mpirecvbuf
allocate(mpirecvbuf(1:lx1,1:lx2,1:lx3))
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit1)
  x%e1(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit2)
  x%e2(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tageunit3)
  x%e3(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tager)
  x%er(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tagetheta)
  x%etheta(:,:,:,icomp)=mpirecvbuf
end do
do icomp=1,3
  call bcast_recv(mpirecvbuf,tagephi)
  x%ephi(:,:,:,icomp)=mpirecvbuf
end do
end block

call bcast_recv(x%r,tagr)
call bcast_recv(x%theta,tagtheta)
call bcast_recv(x%phi,tagphi)


!COUNT THE NUMBER OF NULL GRID POINTS AND GENERATE A LIST OF NULL INDICES FOR LATER USE
x%lnull=0;
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%lnull=x%lnull+1
      end if
    end do
  end do
end do
allocate(x%inull(1:x%lnull,1:3))

icount=1
do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      if(x%nullpts(ix1,ix2,ix3)>0.5d0) then
        x%inull(icount,:)=[ix1,ix2,ix3]
        icount=icount+1
      end if
    end do
  end do
end do
print *, 'Done computing null grid points...  Process:  ',myid,' has:  ',x%lnull


!! COMPUTE DIFFERENTIAL DISTANCES ALONG EACH DIRECTION
block
real(wp), dimension(:,:,:), allocatable :: tmpdx

allocate(x%dl1i(1:lx1,1:lx2,1:lx3),x%dl2i(1:lx1,1:lx2,1:lx3),x%dl3i(1:lx1,1:lx2,1:lx3),tmpdx(1:lx1,1:lx2,1:lx3))

tmpdx=spread(spread(x%dx1i,2,lx2),3,lx3)
x%dl1i=tmpdx*x%h1(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx2i,1,lx1),3,lx3)
x%dl2i=tmpdx*x%h2(1:lx1,1:lx2,1:lx3)
tmpdx=spread(spread(x%dx3i,1,lx1),2,lx2)
x%dl3i=tmpdx*x%h3(1:lx1,1:lx2,1:lx3)
end block

end subroutine read_grid_workers

end submodule grid_read
