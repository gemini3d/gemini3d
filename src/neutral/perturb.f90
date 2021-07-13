submodule (neutral) perturb

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

use grid, only : gridflag
use reader, only : get_neutral2, get_neutral3
use timeutils, only : dateinc, date_filename

use mpimod, only: mpi_realprec, mpi_integer, mpi_comm_world, mpi_status_ignore, &
tag=>gemini_mpi

implicit none (type, external)

external :: mpi_send, mpi_recv

interface ! proj.f90
  module subroutine gridproj_dneu2D(cfg,flagcart,x)
    type(gemini_cfg), intent(in) :: cfg
    logical, intent(in) :: flagcart
    class(curvmesh), intent(inout) :: x
  end subroutine gridproj_dneu2D

  module subroutine gridproj_dneu3D(cfg,x)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(inout) :: x
  end subroutine gridproj_dneu3D
end interface


interface ! interp.f90
  module subroutine spaceinterp_dneu2D(flagcart)
    logical, intent(in) :: flagcart
  end subroutine spaceinterp_dneu2D

  module subroutine spaceinterp_dneu3D()
  end subroutine spaceinterp_dneu3D

  module subroutine timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(inout) :: dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow
    !! intent(out)
  end subroutine timeinterp_dneu
end interface

contains

module procedure neutral_perturb
! subroutine neutral_perturb(cfg,dt,dtneu,t,ymd,UTsec,x,v2grid,v3grid,nn,Tn,vn1,vn2,vn3)

select case (cfg%interptype)
case (0)
!! cartesian interpolation
!! drho inputs (radial distance) interpreted as dy (horizontal distance)
!  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,neudir,drhon,dzn,meanlat,meanlong,x,nn,Tn,vn1,vn2,vn3)
  call neutral_perturb_cart(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
case (1)
!! axisymmetric interpolation
  call neutral_perturb_axisymm(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
case (3)
!! 3D interpolation drhon is takent to be dyn (northward distance)
  call neutral_perturb_3D(dt,dtneu,t,ymd,UTsec,cfg,x,nn,Tn,vn1,vn2,vn3)
case default
  error stop '...Invalid interpolation type specified from input file...'
end select

!Add interpolated perturbations to reference atmosphere arrays
call neutral_update(nn,Tn,vn1,vn2,vn3,v2grid,v3grid)

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

class(curvmesh), intent(inout) :: x
!! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(inout) :: nn
!! intent(out)
!! neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
!! intent(out)

integer :: ix1,ix2,ix3,iid!,irhon,izn
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t+dt/2 >= tnext .or. t < 0) then   !negative time means that we need to load the first frame

  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first time we've tried to load neutral simulation data, should we check for a previous neutral file to load??? or just assume everything starts at zero?  This needs to somehow check for an existing file under certain conditiosn, maybe if it==1???  Actually we don't even need that we can just check that the neutral grid is allocated (or not)
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
  if (mpi_cfg%myid==0 .and. debug) then
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

class(curvmesh), intent(inout) :: x
!! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(inout) :: nn
!! intent(out)
!! neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
!! intent(out)

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t + dt/2 >= tnext .or. t < 0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first time we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !> Create a neutral grid, do some allocations and projections
    call gridproj_dneu2D(cfg,.true.,x)
    !! set true to denote Cartesian...
  end if

  !> Read in neutral data from a file
  call read_dneu2D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp,.true.)

  !> Spatial interpolatin for the frame we just read in
  if (mpi_cfg%myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
  end if
  call spaceinterp_dneu2D(.true.)

  !> UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+dtneu
  UTsecnext=UTsectmp
  ymdnext=ymdtmp

end if

!Interpolation in time
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

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

class(curvmesh), intent(inout) :: x
!! grid structure  (inout because we want to be able to deallocate unit vectors once we are done with them)
real(wp), dimension(:,:,:,:), intent(inout) :: nn
!! intent(out)
!! neutral params interpolated to plasma grid at requested time
real(wp), dimension(:,:,:), intent(inout) :: Tn,vn1,vn2,vn3
!! intent(out)

integer :: ix1,ix2,ix3,iid
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

real(wp) :: starttime,endtime


!CHECK WHETHER WE NEED TO LOAD A NEW FILE
if (t + dt/2 >= tnext .or. t < 0) then
  !IF FIRST LOAD ATTEMPT CREATE A NEUTRAL GRID AND COMPUTE GRID SITES FOR IONOSPHERIC GRID.  Since this needs an input file, I'm leaving it under this condition here
  if (.not. allocated(zn)) then     !means this is the first time we've tried to load neutral simulation data, should we check for a previous neutral file to load???
    !initialize dates
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    !Create a neutral grid, do some allocations and projections
    if (mpi_cfg%myid==0 .and. debug) then
      print*, 'Creating a neutral grid...'
    end if
    call gridproj_dneu3D(cfg,x)
  end if

  !Read in neutral data from a file
  if (mpi_cfg%myid==0 .and. debug) then
    print*, 'Reading in data from neutral file'
    call cpu_time(starttime)
  end if
  call read_dneu3D(tprev,tnext,t,dtneu,dt,cfg%sourcedir,ymdtmp,UTsectmp)
  if (mpi_cfg%myid==0 .and. debug) then
    call cpu_time(endtime)
    print*, 'Neutral data input required time:  ',endtime-starttime
  end if

  !Spatial interpolation for the frame we just read in
  if (mpi_cfg%myid==0 .and. debug) then
    print *, 'Spatial interpolation and rotation of vectors for date:  ',ymdtmp,' ',UTsectmp
    call cpu_time(starttime)
  end if
  call spaceinterp_dneu3D()
  if (mpi_cfg%myid==0 .and. debug) then
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
if (mpi_cfg%myid==0 .and. debug) then
  print*, 'Interpolating in time'
end if
call timeinterp_dneu(t,dt,dNOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow,dTninow)

end subroutine neutral_perturb_3D


subroutine read_dneu2D(tprev,tnext,t,dtneu,dt,neudir,ymdtmp,UTsectmp,flagcart)

! This subroutine reads in neutral frame data, if necessary and puts the results
! in a module-scope variable for later use
! This version is for 2D source neutral data

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt
!! times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir
!! directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp
!! storage space for incrementing date without overwriting ymdnext...
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


if (mpi_cfg%myid==0) then    !root
  !read in the data from file
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
  do iid=1,mpi_cfg%lid-1
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
if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
  print*, 'neutral data size:  ',lhorzn,lzn, mpi_cfg%lid
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

real(wp), intent(in) :: tprev,tnext,t,dtneu,dt
!! times of previous input frame,next frame and then current time
character(*), intent(in) :: neudir
!! directory where neutral simulation data is kept
integer, dimension(3), intent(out) :: ymdtmp
!! storage space for incrementing date without overwriting ymdnext...
real(wp), intent(out) :: UTsectmp


integer :: iid,ierr
integer :: lhorzn                        !number of horizontal grid points
real(wp), dimension(:,:,:), allocatable :: parmtmp    !temporary resizable space for subgrid neutral data


lhorzn=lyn


if (mpi_cfg%myid==0) then    !root
  !read in the data from file
  ymdtmp=ymdnext
  UTsectmp=UTsecnext
  call dateinc(dtneu,ymdtmp,UTsectmp)                !get the date for "next" params

  !FIXME: we probably need to read in and distribute the input parameters one at a time to reduce memory footprint...
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


  ! FIXME: should we loop through the workers for each parameter first? that way we can overwrite that parameter and also not worry about worker 2 having to block until worker 1 receives all of its data...  The issue is that parmtmp changes size for each worker...   Possibly allocated a single linear buffer that can hold the full data and pack that with the data to be sent to avoid repeated allocation and deallocations...  Can we also interleave interpolation with the message passing???
  !in the 3D case we cannot afford to send full grid data and need to instead use neutral subgrid splits defined earlier
  do iid=1,mpi_cfg%lid-1
    allocate(parmtmp(lzn,slabsizes(iid,1),slabsizes(iid,2)))    !get space for the parameters for this worker

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


if (mpi_cfg%myid==mpi_cfg%lid/2 .and. debug) then
  print*, 'neutral data size:  ',mpi_cfg%myid,lzn,lxn,lyn
  print *, 'Min/max values for dnO:  ',mpi_cfg%myid,minval(dnO),maxval(dnO)
  print *, 'Min/max values for dnN:  ',mpi_cfg%myid,minval(dnN2),maxval(dnN2)
  print *, 'Min/max values for dnO2:  ',mpi_cfg%myid,minval(dnO2),maxval(dnO2)
  print *, 'Min/max values for dvnx:  ',mpi_cfg%myid,minval(dvnx),maxval(dvnx)
  print *, 'Min/max values for dvnrho:  ',mpi_cfg%myid,minval(dvnrho),maxval(dvnrho)
  print *, 'Min/max values for dvnz:  ',mpi_cfg%myid,minval(dvnz),maxval(dvnz)
  print *, 'Min/max values for dTn:  ',mpi_cfg%myid,minval(dTn),maxval(dTn)
!  print*, 'coordinate ranges:  ',minval(zn),maxval(zn),minval(rhon),maxval(rhon),minval(zi),maxval(zi),minval(rhoi),maxval(rhoi)
end if

end subroutine read_dneu3D

end submodule perturb
