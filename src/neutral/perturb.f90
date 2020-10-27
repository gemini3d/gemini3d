submodule (neutral) perturb

use grid, only : gridflag
use interpolation, only : interp2, interp3
use mpimod, only : mpi_realprec
use reader, only : get_neutral2, get_neutral3
use timeutils, only : dateinc, date_filename, find_lastdate

use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
myid, lid, tag=>gemini_mpi

implicit none (type, external)

external :: mpi_send, mpi_recv

interface ! proj.f90
module subroutine gridproj_dneu2D(cfg,flagcart,x)
type(gemini_cfg), intent(in) :: cfg
logical, intent(in) :: flagcart
type(curvmesh), intent(inout) :: x
end subroutine gridproj_dneu2D

module subroutine gridproj_dneu3D(cfg,x)
type(gemini_cfg), intent(in) :: cfg
type(curvmesh), intent(inout) :: x
end subroutine gridproj_dneu3D
end interface

contains

module procedure neutral_perturb
! subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)

if (cfg%interptype==0) then                           !cartesian interpolation drho inputs (radial distance) will be interpreted as dy (horizontal distance)
!  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,neudir,drhon,dzn,meanlat,meanlong,x,nn,Tn,vn1,vn2,vn3)
  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else if (cfg%interptype==1) then                      !axisymmetric interpolation
  call neutral_perturb_axisymm(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else if (cfg%interptype==3) then                      !3D interpolation drhon is takent to be dyn (northward distance)
  call neutral_perturb_3D(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
else
  error stop '...Invalid interpolation type specified from input file...'
end if

!> adjust for the grid drift speed; note we assume this gets sets to zero if not a Lagrangian grid
vn2=vn2-v2grid
vn3=vn3-v3grid

end procedure neutral_perturb


!FOR CONSISTENCY I'D LIKE TO STRUCTURE NEUTRAL PERTURB OPERATIONS LIKE GRAVITY IS HANDLED IN THE GRID MODULE, I.E. HAVE AN EXPLICIT CONSTRUCTORS/DESTRUCTOR TYPE ROUTINE THAT HANDLES ALLOCATION AND DEALLOCATION, WHICH WILL CLEAN UP THE NEUTRAL_PERTURB SUBROUTINE, I.E. REMOVE ALLOCATES OF PERSISTENT MODULE VARIABLES.
subroutine neutral_perturb_axisymm(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CYLINDRICAL COORDINATES.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid!,irhon,izn
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2._wp >= tnext .or. t < 0) then   !negative time means that we need to load the first frame

  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load??? or just assume everything starts at zero?  This needs to somehow check for an existing file under certain conditiosn, maybe if it==1???  Actually we don't even need that we can just check that the neutral grid is allocated (or not)
    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    call gridproj_dneu2D(cfg,.false.,x)    !set false to denote not Cartesian; here and below...
  end if

  !Read in neutral data from a file
  call read_dneu2D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp,.false.)

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
  end if
  call spaceinterp_dneu2D(.false.)

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if !done loading frame data...

!Interpolation in time
call timeinterp_dneu(t,dt,dnOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!Add interpolated perturbations to reference atmosphere arrays
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_axisymm


!! THIS SHARES SO MUCH CODE WITH THE AXISYMMETRIC VERSION THAT THEY SHOULD PROBABLY BE COMBINED
subroutine neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER.
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CARTESIAN COORDINATES.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2d0 >= tnext .or. t <= 0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    tprev=t
    tnext=t+cfg%dtneu

    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    call gridproj_dneu2D(cfg,.true.,x)    !set true to denote Cartesian...
  end if

  !Read in neutral data from a file
  call read_dneu2D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp,.true.)

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
  end if
  call spaceinterp_dneu2D(.true.)

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if

!Interpolation in time
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!> update neutral atmosphere
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_cart


subroutine neutral_perturb_3D(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)

!------------------------------------------------------------
!-------COMPUTE NEUTRAL PERTURBATIONS FOR THIS TIME STEP.  ADD
!-------THEM TO MSIS PERTURBATIONS TO GET ABSOLUTE VALUES FOR
!-------EACH PARAMETER.
!-------
!-------THIS VERSION ASSUMES THE INPUT NEUTRAL DATA ARE IN
!-------CARTESIAN COORDINATES AND IN THREE DIMENSIONS.
!------------------------------------------------------------

real(wp), intent(in) :: dt,dtneu
real(wp), intent(in) :: t
integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
type(gemini_cfg), intent(in) :: cfg

type(curvmesh), intent(inout) :: x         !grid structure  (inout becuase we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(out) :: nn   !neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(out) :: Tn,vn1,vn2,vn3

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

real(wp) :: starttime,endtime


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t + dt/2._wp >= tnext .or. t <= 0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first tiem we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    tprev=t
    tnext=t+cfg%dtneu

    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    if (myid==0 .and. debug) then
      print*, 'Creating a neutral grid...'
    end if
    call gridproj_dneu3D(cfg,x)
  end if

  !Read in neutral data from a file
  if (myid==0 .and. debug) then
    print*, 'Reading in data from neutral file'
    call cpu_time(starttime)
  end if
  call read_dneu3D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp)
  if (myid==0 .and. debug) then
    call cpu_time(endtime)
    print*, 'Neutral data input required time:  ',endtime-starttime
  end if

  !Spatial interpolatin for the frame we just read in
  if (myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
    call cpu_time(starttime)
  end if
  call spaceinterp_dneu3D()
  if (myid==0 .and. debug) then
    call cpu_time(endtime)
    print*, 'Spatial interpolation in 3D took time:  ',endtime-starttime
  end if

  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if

!Interpolation in time
if (myid==0 .and. debug) then
  print*, 'Interpolating in time'
end if
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!> update neutral atmosphere
call neutral_update(nn,Tn,vn1,vn2,vn3)

end subroutine neutral_perturb_3D


subroutine read_dneu2D(tprev,tnext,t,dtneu,dt,neudir,ymdtmp,UTsectmp,flagcart)

! This subroutine reads in neutral frame data, if necessary and puts the results
! in a module-scope variable for later use
! This version is for 2D source neutral data

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt        !times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir                    !directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp          !storage space for incrementing date without overwriting ymdnext...
real(wp), intent(out) :: UTsectmp
logical, intent(in) :: flagcart

integer :: iid,ierr
integer :: lhorzn
!! number of horizontal grid points


if (flagcart) then
  lhorzn=lyn
else
  lhorzn=lrhon
end if


if (myid==0) then    !root
  !read in the data from file
  if(debug) print *, 'neutral.f90:read_dneu2D: tprev,tnow,tnext:  ',tprev,t+dt/2,tnext
  ymdtmp=ymdnext
  UTsectmp=UTsecnext
  call dateinc(dtneu,ymdtmp,UTsectmp)    !get the date for "next" params

  call get_neutral2(date_filename(neudir,ymdtmp,UTsectmp), &
    dnO,dnN2,dnO2,dvnrho,dvnz,dTn)

  if (debug) then
    print *, 'Min/max values for dnO:  ',minval(dnO),maxval(dnO)
    print *, 'Min/max values for dnN:  ',minval(dnN2),maxval(dnN2)
    print *, 'Min/max values for dnO:  ',minval(dnO2),maxval(dnO2)
    print *, 'Min/max values for dvnrho:  ',minval(dvnrho),maxval(dvnrho)
    print *, 'Min/max values for dvnz:  ',minval(dvnz),maxval(dvnz)
    print *, 'Min/max values for dTn:  ',minval(dTn),maxval(dTn)
  endif

  if (.not. all(ieee_is_finite(dnO))) error stop 'dnO: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnN2))) error stop 'dnN2: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnO2))) error stop 'dnO2: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnrho))) error stop 'dvnrho: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnz))) error stop 'dvnz: non-finite value(s)'
  if (.not. all(ieee_is_finite(dTn))) error stop 'dTn: non-finite value(s)'

  !send a full copy of the data to all of the workers
  do iid=1,lid-1
    call mpi_send(dnO,lhorzn*lzn,mpi_realprec,iid,tag%dnO,MPI_COMM_WORLD,ierr)
    call mpi_send(dnN2,lhorzn*lzn,mpi_realprec,iid,tag%dnN2,MPI_COMM_WORLD,ierr)
    call mpi_send(dnO2,lhorzn*lzn,mpi_realprec,iid,tag%dnO2,MPI_COMM_WORLD,ierr)
    call mpi_send(dTn,lhorzn*lzn,mpi_realprec,iid,tag%dTn,MPI_COMM_WORLD,ierr)
    call mpi_send(dvnrho,lhorzn*lzn,mpi_realprec,iid,tag%dvnrho,MPI_COMM_WORLD,ierr)
    call mpi_send(dvnz,lhorzn*lzn,mpi_realprec,iid,tag%dvnz,MPI_COMM_WORLD,ierr)
  end do
else     !workers
  !receive a full copy of the data from root
  call mpi_recv(dnO,lhorzn*lzn,mpi_realprec,0,tag%dnO,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnN2,lhorzn*lzn,mpi_realprec,0,tag%dnN2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnO2,lhorzn*lzn,mpi_realprec,0,tag%dnO2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dTn,lhorzn*lzn,mpi_realprec,0,tag%dTn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnrho,lhorzn*lzn,mpi_realprec,0,tag%dvnrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnz,lhorzn*lzn,mpi_realprec,0,tag%dvnz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


!DO SPATIAL INTERPOLATION OF EACH PARAMETER (COULD CONSERVE SOME MEMORY BY NOT STORING DVNRHOIPREV AND DVNRHOINEXT, ETC.)
if (myid==lid/2 .and. debug) then
  print*, 'neutral data size:  ',lhorzn,lzn,lid
  print *, 'Min/max values for dnO:  ',minval(dnO),maxval(dnO)
  print *, 'Min/max values for dnN:  ',minval(dnN2),maxval(dnN2)
  print *, 'Min/max values for dnO:  ',minval(dnO2),maxval(dnO2)
  print *, 'Min/max values for dvnrho:  ',minval(dvnrho),maxval(dvnrho)
  print *, 'Min/max values for dvnz:  ',minval(dvnz),maxval(dvnz)
  print *, 'Min/max values for dTn:  ',minval(dTn),maxval(dTn)
  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
end if

end subroutine read_dneu2D


subroutine read_dneu3D(tprev,tnext,t,dtneu,dt,neudir,ymdtmp,UTsectmp)

! This subroutine reads in neutral frame data, if necessary and puts the results
! in a module-scope variable for later use
! this version is for 3D neutral source data

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt        !times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir                    !directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp          !storage space for incrementing date without overwriting ymdnext...
real(wp), intent(out) :: UTsectmp


integer :: iid,ierr
integer :: lhorzn                        !number of horizontal grid points
real(wp), dimension(:,:,:), allocatable :: parmtmp    !temporary resizable space for subgrid neutral data


lhorzn=lyn


if (myid==0) then    !root
  !read in the data from file
  if(debug) print *, 'tprev,tnow,tnext:  ',tprev,t+dt/2d0,tnext
  ymdtmp=ymdnext
  UTsectmp=UTsecnext
  call dateinc(dtneu,ymdtmp,UTsectmp)                !get the date for "next" params

  call get_neutral3(date_filename(neudir,ymdtmp,UTsectmp), &
    dnOall,dnN2all,dnO2all,dvnxall,dvnrhoall,dvnzall,dTnall)

  if (debug) then
    print *, 'Min/max values for dnOall:  ',minval(dnOall),maxval(dnOall)
    print *, 'Min/max values for dnNall:  ',minval(dnN2all),maxval(dnN2all)
    print *, 'Min/max values for dnO2all:  ',minval(dnO2all),maxval(dnO2all)
    print *, 'Min/max values for dvnxall:  ',minval(dvnxall),maxval(dvnxall)
    print *, 'Min/max values for dvnrhoall:  ',minval(dvnrhoall),maxval(dvnrhoall)
    print *, 'Min/max values for dvnzall:  ',minval(dvnzall),maxval(dvnzall)
    print *, 'Min/max values for dTnall:  ',minval(dTnall),maxval(dTnall)
  endif

  if (.not. all(ieee_is_finite(dnOall))) error stop 'dnOall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnN2all))) error stop 'dnN2all: non-finite value(s)'
  if (.not. all(ieee_is_finite(dnO2all))) error stop 'dnO2all: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnxall))) error stop 'dvnxall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnrhoall))) error stop 'dvnrhoall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dvnzall))) error stop 'dvnzall: non-finite value(s)'
  if (.not. all(ieee_is_finite(dTnall))) error stop 'dTnall: non-finite value(s)'


  !in the 3D case we cannnot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
  do iid=1,lid-1
    allocate(parmtmp(lzn,slabsizes(iid,1),slabsizes(iid,2)))    !get space for the parameter for this worker

    parmtmp=dnOall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnO,MPI_COMM_WORLD,ierr)

    parmtmp=dnN2all(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnN2,MPI_COMM_WORLD,ierr)

    parmtmp=dnO2all(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dnO2,MPI_COMM_WORLD,ierr)

    parmtmp=dTnall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dTn,MPI_COMM_WORLD,ierr)

    parmtmp=dvnrhoall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnrho,MPI_COMM_WORLD,ierr)

    parmtmp=dvnzall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnz,MPI_COMM_WORLD,ierr)

    parmtmp=dvnxall(1:lzn,indx(iid,3):indx(iid,4),indx(iid,5):indx(iid,6))
    call mpi_send(parmtmp,lzn*slabsizes(iid,1)*slabsizes(iid,2),mpi_realprec,iid,tag%dvnx,MPI_COMM_WORLD,ierr)

    deallocate(parmtmp)
  end do

  !root needs to grab its piece of data
  dnO=dnOall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dnN2=dnN2all(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dnO2=dnO2all(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dTn=dTnall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnrho=dvnrhoall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnz=dvnzall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
  dvnx=dvnxall(1:lzn,indx(0,3):indx(0,4),indx(0,5):indx(0,6))
else     !workers
  !receive a subgrid copy of the data from root
  call mpi_recv(dnO,lzn*lxn*lyn,mpi_realprec,0,tag%dnO,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnN2,lzn*lxn*lyn,mpi_realprec,0,tag%dnN2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dnO2,lzn*lxn*lyn,mpi_realprec,0,tag%dnO2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dTn,lzn*lxn*lyn,mpi_realprec,0,tag%dTn,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnrho,lzn*lxn*lyn,mpi_realprec,0,tag%dvnrho,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnz,lzn*lxn*lyn,mpi_realprec,0,tag%dvnz,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call mpi_recv(dvnx,lzn*lxn*lyn,mpi_realprec,0,tag%dvnx,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if


if (myid==lid/2 .and. debug) then
  print*, 'neutral data size:  ',myid,lzn,lxn,lyn
  print *, 'Min/max values for dnO:  ',myid,minval(dnO),maxval(dnO)
  print *, 'Min/max values for dnN:  ',myid,minval(dnN2),maxval(dnN2)
  print *, 'Min/max values for dnO2:  ',myid,minval(dnO2),maxval(dnO2)
  print *, 'Min/max values for dvnx:  ',myid,minval(dvnx),maxval(dvnx)
  print *, 'Min/max values for dvnrho:  ',myid,minval(dvnrho),maxval(dvnrho)
  print *, 'Min/max values for dvnz:  ',myid,minval(dvnz),maxval(dvnz)
  print *, 'Min/max values for dTn:  ',myid,minval(dTn),maxval(dTn)
!  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
end if

end subroutine read_dneu3D

subroutine spaceinterp_dneu2D(flagcart)

!must take into account the type of interpolation that is being done
logical, intent(in) :: flagcart

integer :: lhorzn
real(wp), dimension(:,:), allocatable :: tmpinterp
real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data


!Array for packing neutral data
if (flagcart) then
  lhorzn=lyn
else
  lhorzn=lrhon
end if
allocate(tmpinterp(lzn,lhorzn))


if(flagcart) then
  tmpinterp=dnO(:,1,:)                    !pack into 2D array for interp2
  parami=interp2(zn,yn,tmpinterp,zi,yi)         !interp to temp var.
  dnOiprev=dnOinext                       !save new previous
  dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

  tmpinterp=dnN2(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dnN2iprev=dnN2inext
  dnN2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dnO2(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dnO2iprev=dnO2inext
  dnO2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnrho(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
  dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnz(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dvnziprev=dvnzinext
  dvnzinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dTn(:,1,:)
  parami=interp2(zn,yn,tmpinterp,zi,yi)
  dTniprev=dTninext
  dTninext=reshape(parami,[lx1,lx2,lx3])
else
  tmpinterp=dnO(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)     !interp to temp var.
  dnOiprev=dnOinext                       !save new previous
  dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

  tmpinterp=dnN2(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dnN2iprev=dnN2inext
  dnN2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dnO2(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dnO2iprev=dnO2inext
  dnO2inext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnrho(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dvnrhoiprev=dvnrhoinext
  dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dvnz(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dvnziprev=dvnzinext
  dvnzinext=reshape(parami,[lx1,lx2,lx3])

  tmpinterp=dTn(:,1,:)
  parami=interp2(zn,rhon,tmpinterp,zi,rhoi)
  dTniprev=dTninext
  dTninext=reshape(parami,[lx1,lx2,lx3])
end if

!MORE DIAG
if (myid==lid/2) then
  print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvrhoi:  ',minval(dvnrhoinext),maxval(dvnrhoinext)
  print *, 'Min/max values for dvnzi:  ',minval(dvnzinext),maxval(dvnzinext)
  print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
end if


!ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
!structure)
dvn1iprev=dvn1inext   !save the old data
dvn2iprev=dvn2inext
dvn3iprev=dvn3inext
if(flagcart) then
  dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
  dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2
  dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3
else
  dvn1inext=dvnrhoinext*proj_erhop_e1+dvnzinext*proj_ezp_e1    !apply projection to complete rotation into dipole coordinates
  dvn2inext=dvnrhoinext*proj_erhop_e2+dvnzinext*proj_ezp_e2
  dvn3inext=dvnrhoinext*proj_erhop_e3+dvnzinext*proj_ezp_e3
end if

!MORE DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dnOi:  ',minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvn1i:  ',minval(dvn1inext),maxval(dvn1inext)
  print *, 'Min/max values for dvn2i:  ',minval(dvn2inext),maxval(dvn2inext)
  print *, 'Min/max values for dvn3i:  ',minval(dvn3inext),maxval(dvn3inext)
  print *, 'Min/max values for dTni:  ',minval(dTninext),maxval(dTninext)
end if


!CLEAR ALLOCATED VARS
deallocate(tmpinterp)

end subroutine spaceinterp_dneu2D


subroutine spaceinterp_dneu3D()

!performs spatial interpolation for 3D input neutral data from MAGIC or some other source

real(wp), dimension(lx1*lx2*lx3) :: parami    !work array for temp storage of interpolated data, note sizes taken from grid module data


!INTERPOLATE IN THREE DIMENSIONS
parami=interp3(zn,xn,yn,dnO,zi,xi,yi)         !interp to temp var.
dnOiprev=dnOinext                       !save new previous
dnOinext=reshape(parami,[lx1,lx2,lx3])  !overwrite next with new interpolated input

parami=interp3(zn,xn,yn,dnN2,zi,xi,yi)
dnN2iprev=dnN2inext
dnN2inext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dnO2,zi,xi,yi)
dnO2iprev=dnO2inext
dnO2inext=reshape(parami,[lx1,lx2,lx3])

!ZZZ - do we want to make dvnrho-->dvny???
parami=interp3(zn,xn,yn,dvnrho,zi,xi,yi)
dvnrhoiprev=dvnrhoinext    !interpreted as y-component in this (cartesian) function
dvnrhoinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnz,zi,xi,yi)
dvnziprev=dvnzinext
dvnzinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dvnx,zi,xi,yi)
dvnxiprev=dvnxinext
dvnxinext=reshape(parami,[lx1,lx2,lx3])

parami=interp3(zn,xn,yn,dTn,zi,xi,yi)
dTniprev=dTninext
dTninext=reshape(parami,[lx1,lx2,lx3])


!MORE DIAG
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dnOi:  ',myid,minval(dnOinext),maxval(dnOinext)
  print *, 'Min/max values for dnN2i:  ',myid,minval(dnN2inext),maxval(dnN2inext)
  print *, 'Min/max values for dnO2i:  ',myid,minval(dnO2inext),maxval(dnO2inext)
  print *, 'Min/max values for dvrhoi:  ',myid,minval(dvnrhoinext),maxval(dvnrhoinext)
  print *, 'Min/max values for dvnzi:  ',myid,minval(dvnzinext),maxval(dvnzinext)
  print *, 'Min/max values for dvnxi:  ',myid,minval(dvnxinext),maxval(dvnxinext)
  print *, 'Min/max values for dTni:  ',myid,minval(dTninext),maxval(dTninext)
end if


!ROTATE VECTORS INTO X1 X2 DIRECTIONS (Need to include unit vectors with grid
!structure)
dvn1iprev=dvn1inext   !save the old data for the rotated vectors
dvn2iprev=dvn2inext
dvn3iprev=dvn3inext
dvn1inext=dvnrhoinext*proj_eyp_e1+dvnzinext*proj_ezp_e1+dvnxinext*proj_exp_e1    !apply projection to complete rotation into dipole coordinates; drhoi interpreted here at teh y component (northward)
dvn2inext=dvnrhoinext*proj_eyp_e2+dvnzinext*proj_ezp_e2+dvnxinext*proj_exp_e2
dvn3inext=dvnrhoinext*proj_eyp_e3+dvnzinext*proj_ezp_e3+dvnxinext*proj_exp_e3


!MORE DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'Min/max values for dvn1i:  ',myid,minval(dvn1inext),maxval(dvn1inext)
  print *, 'Min/max values for dvn2i:  ',myid,minval(dvn2inext),maxval(dvn2inext)
  print *, 'Min/max values for dvn3i:  ',myid,minval(dvn3inext),maxval(dvn3inext)
end if

end subroutine spaceinterp_dneu3D


subroutine timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

!interpolation in time - no sensitive to dimensionality of the input neutral data so this can be
! the same for 2D vs. 3D

real(wp), intent(in) :: t,dt
real(wp), dimension(:,:,:), intent(out) :: dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow

integer :: ix1,ix2,ix3
real(wp) :: slope

do ix3=1,lx3
  do ix2=1,lx2
    do ix1=1,lx1
      slope=(dnOinext(ix1,ix2,ix3)-dnOiprev(ix1,ix2,ix3))/(tnext-tprev)
      dnOinow(ix1,ix2,ix3)=dnOiprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnN2inext(ix1,ix2,ix3)-dnN2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnN2inow(ix1,ix2,ix3)=dnN2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dnO2inext(ix1,ix2,ix3)-dnO2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dnO2inow(ix1,ix2,ix3)=dnO2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn1inext(ix1,ix2,ix3)-dvn1iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn1inow(ix1,ix2,ix3)=dvn1iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn2inext(ix1,ix2,ix3)-dvn2iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn2inow(ix1,ix2,ix3)=dvn2iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dvn3inext(ix1,ix2,ix3)-dvn3iprev(ix1,ix2,ix3))/(tnext-tprev)
      dvn3inow(ix1,ix2,ix3)=dvn3iprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)

      slope=(dTninext(ix1,ix2,ix3)-dTniprev(ix1,ix2,ix3))/(tnext-tprev)
      dTninow(ix1,ix2,ix3)=dTniprev(ix1,ix2,ix3)+slope*(t+dt/2.0_wp-tprev)
    end do
  end do
end do


!SOME BASIC DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'tprev,t,tnext:  ',myid,tprev,t+dt/2d0,tnext
  print *, 'Min/max values for dnOinow:  ',myid,minval(dnOinow),maxval(dnOinow)
  print *, 'Min/max values for dnN2inow:  ',myid,minval(dnN2inow),maxval(dnN2inow)
  print *, 'Min/max values for dnO2inow:  ',myid,minval(dnO2inow),maxval(dnO2inow)
  print *, 'Min/max values for dvn1inow:  ',myid,minval(dvn1inow),maxval(dvn1inow)
  print *, 'Min/max values for dvn2inow:  ',myid,minval(dvn2inow),maxval(dvn2inow)
  print *, 'Min/max values for dvn3inow:  ',myid,minval(dvn3inow),maxval(dvn3inow)
  print *, 'Min/max values for dTninow:  ',myid,minval(dTninow),maxval(dTninow)
end if

end subroutine timeinterp_dneu


subroutine slabrange(maxzn,ximat,yimat,zimat,sourcemlat,xnrange,ynrange)

!takes in a subgrid and the max altitude of interest for neutral interpolation and then computes
!what the maximum xn and yn will be for that slab
! ZZZ - also this is specific to dipole grids right now...

real(wp), intent(in) :: maxzn
real(wp), dimension(:,:,:), intent(in) :: ximat,yimat,zimat
real(wp), intent(in) :: sourcemlat
real(wp), dimension(2), intent(out) :: xnrange,ynrange     !for min and max

real(wp), dimension(:,:,:), allocatable :: xitmp,yitmp,zitmp
integer :: lx1tmp
integer, dimension(size(ximat,2),size(ximat,3)) :: ix1stmp
integer :: ix1tmp
logical :: flagSH
integer :: ix1


!in what hemisphere is our source?
if (sourcemlat<=0) then
  flagSH=.true.
else
  flagSH=.false.
end if


!peel the grid in half (source hemisphere if closed dipole)
if (gridflag==0) then    !closed dipole grid

  ix1=maxloc(pack(zimat(:,1,1),.true.),1)    !apex is by definition the highest altitude along a given field line
  if (flagSH) then
    lx1tmp=ix1                  !first piece of arrays
  else
    lx1tmp=lx1-ix1    !second (end) piece of arrays
  end if
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers???

  if(flagSH) then    !southern hemisphere
    xitmp=ximat(1:ix1,1:lx2,1:lx3)          !select beginning of the array - the southern half
    yitmp=yimat(1:ix1,1:lx2,1:lx3)
    zitmp=zimat(1:ix1,1:lx2,1:lx3)
  else               !northern hemisphere
      xitmp=ximat(ix1+1:lx1,1:lx2,1:lx3)    !select end half of the array
      yitmp=yimat(ix1+1:lx1,1:lx2,1:lx3)
      zitmp=zimat(ix1+1:lx1,1:lx2,1:lx3)
  end if
else     !this is not an interhemispheric grid so our approach is to just use all of the data
  lx1tmp=lx1
  allocate(xitmp(lx1tmp,lx2,lx3),yitmp(lx1tmp,lx2,lx3),zitmp(lx1tmp,lx2,lx3))   !could this be done more less wastefully with pointers?
  xitmp=ximat(1:lx1,1:lx2,1:lx3)
  yitmp=yimat(1:lx1,1:lx2,1:lx3)
  zitmp=zimat(1:lx1,1:lx2,1:lx3)
!  flagSH=.true.    !treat is as southern, doesn't really matter in this case...
end if


!the min and max x are simply determined by longitude...
xnrange(1)=minval(xitmp)
xnrange(2)=maxval(xitmp)


!situation is more complicated for latitude due to dipole grid, need to determine by L-shell
if (flagSH) then
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0, just use the first longitude slice since they will all have the same L-shell-field line relations
  ynrange(2)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(1)=yitmp(ix1,lx2,1)
else    !things are swapped around in NH
  ix1=minloc(zitmp(:,1,1)-maxzn,1,zitmp(:,1,1)-maxzn>0._wp)    !find the min distance from maxzn subject to constraint that it is > 0; this is the southernmost edge of the neutral slab we need
  ynrange(1)=yitmp(ix1,1,1)
  ix1=minloc(zitmp(:,lx2,1),1,zitmp(:,lx2,1)<0._wp)
  ynrange(2)=yitmp(ix1,lx2,1)
end if

deallocate(xitmp,yitmp,zitmp)

end subroutine slabrange


subroutine range2inds(ranges,zn,xnall,ynall,indices)

!determine where the slab described by ranges falls within the global neutral grid

real(wp), dimension(6), intent(in) :: ranges
real(wp), dimension(:), intent(in) :: zn,xnall,ynall
integer, dimension(6), intent(out) :: indices

real(wp) :: minzn,maxzn,minxn,maxxn,minyn,maxyn
integer :: ixn,iyn


!for clarity
minzn=ranges(1)
maxzn=ranges(2)
minxn=ranges(3)
maxxn=ranges(4)
minyn=ranges(5)
maxyn=ranges(6)


!always use the full z-range
indices(1)=1
indices(2)=lzn


!x-range
ixn=1
do while (ixn<lxnall .and. xnall(ixn)<minxn)
  ixn=ixn+1
end do
indices(3)=max(ixn-1,1)    !just to be sure go back one index so that we cover the min range, don't let index fall below zero
do while (ixn<lxnall .and. xnall(ixn)<maxxn)
  ixn=ixn+1
end do
indices(4)=ixn


!y-range
iyn=1
do while (iyn<lynall .and. ynall(iyn)<minyn)
  iyn=iyn+1
end do
indices(5)=max(iyn-1,1)    !just to be sure go back one index so that we cover the min range
do while (iyn<lynall .and. ynall(iyn)<maxyn)
  iyn=iyn+1
end do
indices(6)=iyn

print*, '!!!!!!!!!!!!!!!!!'
print*, myid
print*, ranges
print*, indices
print*, lxnall,lynall
print*, xnall(indices(3)),xnall(indices(4))
print*, ynall(indices(5)),ynall(indices(6))
print*, '!!!!!!!!!!!!!!!!!'


!corner cases - range is not at all within the neutral grid...  Manifests as both indices being either 1 or lxi, interpolation should zero these out...

end subroutine range2inds



end submodule perturb
