module precipBCs_mod
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use reader, only: get_simsize2, get_precip, get_grid2
use phys_consts, only: pi,wp, debug
use grid, only : lx1,lx2,lx3,lx3all
use mesh, only: curvmesh
use interpolation, only : interp1,interp2
use timeutils, only : dateinc, date_filename, find_lastdate
use mpimod, only: mpi_integer, mpi_comm_world, mpi_status_ignore, &
lid, mpi_realprec, myid, tag=>gemini_mpi
use config, only: gemini_cfg

implicit none (type, external)
private
public :: clear_precip_fileinput, precipBCs_fileinput, precipBCs, init_precipinput

external :: mpi_send, mpi_recv

!ALL OF THE FOLLOWING MODULE-SCOPE ARRAYS ARE USED FOR INTERPOLATING PRECIPITATION INPUT FILES (IF USED)
real(wp), dimension(:), allocatable, private :: mlonp
real(wp), dimension(:), allocatable, private :: mlatp    !coordinates of precipitation data
integer, private :: llon,llat

real(wp), dimension(:,:), allocatable, private :: Qp,E0p    !total energy flux and char. energy of input data
real(wp), dimension(:), allocatable, private :: precdatp    !needed when a 1D interpolation is to be done, i.e. when there is 1D sourde data

real(wp), dimension(:), allocatable, private :: mloni    !flat list of mlat,mlon locations on grid that we need to interpolate onto
real(wp), dimension(:), allocatable, private :: mlati

real(wp), dimension(:,:), allocatable, private :: Qiprev,E0iprev,Qinext,E0inext    !interpolated in space energy flux and char. en.

integer, dimension(3), private :: ymdprev,ymdnext   !dates for interpolated data
real(wp), private :: UTsecprev,UTsecnext
real(wp), private :: tprev=0._wp,tnext=0._wp


contains


subroutine init_precipinput(dt,t,cfg,ymd,UTsec,x)

!> initialize variables to hold input file precipitation information, must be called by all workers at the same time

real(wp), intent(in) :: dt,t
type(gemini_cfg), intent(in) :: cfg
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
type(curvmesh), intent(in) :: x

real(wp), dimension(1:x%lx2,1:x%lx3,2) :: W0,PhiWmWm2    ! these are only worker-sized, hardcoded 2 precipitation populations...
integer, dimension(3) :: ymdtmp
real(wp) :: UTsectmp

if (cfg%flagprecfile==1) then    !all workers must have this info
  !! find the last input data preceding the milestone/initial condition that we start with
  call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,cfg%dtE0,ymdtmp,UTsectmp)

  if (myid==0) print*, '!!!Attmpting to prime precipitation input files...',ymdtmp,UTsectmp
  !! back up by two dtprec to so that when we run the fileinput twice we end up with tprev corresponding
  !   to the first time step
  tprev=UTsectmp-UTsec-2._wp*cfg%dtprec
  tnext=tprev+cfg%dtprec
  call precipBCs_fileinput(dt,tnext+cfg%dtprec/2._wp,cfg,ymdtmp,UTsectmp-cfg%dtprec,x,W0,PhiWmWm2)
  !time that we interpolate to doesn't matter but it needs to trigger a new read...

  if (myid==0) print*, 'Now loading initial next file for precipitation input...'
  !! now shift next->prev and load a new one corresponding to the first simulation time step
  call precipBCs_fileinput(dt,0._wp,cfg,ymdtmp,UTsectmp,x,W0,PhiWmWm2)
end if

end subroutine init_precipinput


!subroutine set_tmarkers(tprevtmp,tnexttmp)
!
!!> set tprev and tnext to a specified value
!
!tprev=tprevtmp
!tnext=tnexttmp
!
!end subroutine


subroutine precipBCs_fileinput(dt,t,cfg,ymd,UTsec,x,W0,PhiWmWm2)

real(wp), intent(in) :: dt
real(wp), intent(in) :: t
type(gemini_cfg), intent(in) :: cfg
integer, dimension(3), intent(in) :: ymd
!! date for which we wish to calculate perturbations
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: W0,PhiWmWm2
!! last dimension is the number of particle populations
type(curvmesh), intent(in) :: x

integer :: ios, ierr
integer :: iid,iflat,ix2,ix3

real(wp) :: UTsectmp
integer, dimension(3) :: ymdtmp

real(wp), dimension(lx2*lx3) :: parami
real(wp), dimension(lx2,lx3) :: slope,Qinow,E0inow
real(wp) :: W0pk,PhiWpk

UTsectmp = 0

if(t+dt / 2._wp>=tnext .or. t<0._wp) then    !need to load a new file
  if ( .not. allocated(mlonp)) then    !need to read in the grid data from input file
    ymdprev=ymd
    UTsecprev=UTsec
    ymdnext=ymdprev
    UTsecnext=UTsecprev

    if (myid==0) then    !root process
      !READ IN THE GRID
      print '(/,A,/,A)', 'Precipitation input:','--------------------'
      print '(A)', 'READ precipitation size from: ' // cfg%precdir

      call get_simsize2(cfg%precdir, llon=llon, llat=llat)

      print '(A,2I6)', 'Precipitation size: llon,llat:  ',llon,llat
      if (llon < 1 .or. llat < 1) then
        write(stderr,*) 'ERROR: reading ' // cfg%precdir
        error stop 'precipBCs_mod: precipitation grid size must be strictly positive'
      endif

      !MESSAGE WORKERS WITH GRID INFO
      ierr=0
      do iid=1,lid-1
        call mpi_send(llon,1,MPI_INTEGER,iid,tag%llon,MPI_COMM_WORLD,ierr)
        call mpi_send(llat,1,MPI_INTEGER,iid,tag%llat,MPI_COMM_WORLD,ierr)
      end do
      allocate(mlonp(llon),mlatp(llat))    !bit of code duplication with worker code block below...

      if (ierr /= 0) error stop 'mpi_send failed to send grid info'
      !IF WE HAVE SINGLETON DIMENSION THEN ALLOCATE SOME SPACE FOR A TEMP
      !ARRAY FOR INPUTTING INTO INTERP1
      if (llon==1) then
        allocate(precdatp(llat))
      elseif (llat==1) then
        allocate(precdatp(llon))
      end if


      !NOW READ THE GRID
      call get_grid2(cfg%precdir, mlonp, mlatp)

      print '(A,4F9.3)', 'Precipitation mlon,mlat extent:  ',minval(mlonp(:)),maxval(mlonp(:)),minval(mlatp(:)), &
                                                              maxval(mlatp(:))
      if(.not. all(ieee_is_finite(mlonp))) error stop 'mlon must be finite'
      if(.not. all(ieee_is_finite(mlatp))) error stop 'mlat must be finite'

      !NOW SEND THE GRID DATA
      do iid=1,lid-1
        call mpi_send(mlonp,llon,mpi_realprec,iid,tag%mlon,MPI_COMM_WORLD,ierr)
        call mpi_send(mlatp,llat,mpi_realprec,iid,tag%mlat,MPI_COMM_WORLD,ierr)
      end do
    else    !workers
      call mpi_recv(llon,1,MPI_INTEGER,0,tag%llon,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call mpi_recv(llat,1,MPI_INTEGER,0,tag%llat,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      allocate(mlonp(llon),mlatp(llat))

      call mpi_recv(mlonp,llon,mpi_realprec,0,tag%mlon,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call mpi_recv(mlatp,llat,mpi_realprec,0,tag%mlat,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end if

    !SPACE TO STORE INPUT DATA
    allocate(Qp(llon,llat),E0p(llon,llat))
    allocate(Qiprev(lx2,lx3),E0iprev(lx2,lx3),Qinext(lx2,lx3),E0inext(lx2,lx3))
    Qiprev = 0
    E0iprev = 100._wp
    Qinext = 0
    E0inext = 100._wp
    !! these need to be initialized so that something sensible happens at the beginning

    !ALL PROCESSES NEED TO DEFINED THE OPINTS THAT THEY WILL BE INTERPOLATING ONTO
    allocate(mloni(lx2*lx3),mlati(lx2*lx3))
    do ix3=1,lx3
      do ix2=1,lx2
        iflat=(ix3-1)*lx2+ix2
        mlati(iflat)=90._wp-x%theta(lx1,ix2,ix3)*180._wp/pi
        mloni(iflat)=x%phi(lx1,ix2,ix3)*180._wp/pi
      end do
    end do

  end if


  !GRID INFORMATION EXISTS AT THIS POINT SO START READING IN PRECIP DATA
  if (myid==0) then    !only root reads file data
    !read in the data from file
    if(debug) print *, 'precipBCs_mod.f90:precipBCs_fileinput:tprev,tnow,tnext:  ',tprev,t+dt / 2._wp,tnext
    ymdtmp=ymdnext
    UTsectmp=UTsecnext
    call dateinc(cfg%dtprec,ymdtmp,UTsectmp)    !get the date for "next" params

    call get_precip(date_filename(cfg%precdir,ymdtmp,UTsectmp), Qp, E0p)

    if (debug) print *, 'Min/max values for Qp:  ',minval(Qp),maxval(Qp)
    if (debug) print *, 'Min/max values for E0p:  ',minval(E0p),maxval(E0p)

    if(.not. all(ieee_is_finite(Qp))) error stop 'precipBCs_mod.f90:precipBCs_fileinput: Qp must be finite'
    if(any(Qp < 0)) error stop 'precipBCs_mod.f90:precipBCs_fileinput: Qp must be non-negative'
    if(.not. all(ieee_is_finite(E0p))) error stop 'precipBCs_mod.f90:precipBCs_fileinput: E0p must be finite'
    if(any(E0p < 0)) error stop 'precipBCs_mod.f90:precipBCs_fileinput: E0p must be non-negative'

    !send a full copy of the data to all of the workers
    do iid=1,lid-1
      call mpi_send(Qp,llon*llat,mpi_realprec,iid,tag%Qp,MPI_COMM_WORLD,ierr)
      call mpi_send(E0p,llon*llat,mpi_realprec,iid,tag%E0p,MPI_COMM_WORLD,ierr)
    end do
  else     !workers receive data from root
    call mpi_recv(Qp,llon*llat,mpi_realprec,0,tag%Qp,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call mpi_recv(E0p,llon*llat,mpi_realprec,0,tag%E0p,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  end if


  !ALL WORKERS DO SPATIAL INTERPOLATION
  if (myid==0) then
    if (debug) print *, 'Initiating precipitation spatial interpolations for date:  ',ymdtmp,' ',UTsectmp
  end if
  if (llon==1) then    !source data has singleton size in the longitude dimension
    precdatp=Qp(1,:)
    parami=interp1(mlatp,precdatp,mlati)   !will work even for 2D grids, just repeats the data in the lon direction
    Qiprev=Qinext
    Qinext=reshape(parami,[lx2,lx3])

    precdatp=E0p(1,:)
    parami=interp1(mlatp,precdatp,mlati)   !will work even for 2D grids, just repeats the data in the lon direction
    E0iprev=E0inext
    E0inext=reshape(parami,[lx2,lx3])
  elseif (llat==1) then     !source is singleton in lat.
    precdatp=Qp(:,1)
    parami=interp1(mlonp,precdatp,mloni)
    Qiprev=Qinext
    Qinext=reshape(parami,[lx2,lx3])

    precdatp=E0p(:,1)
    parami=interp1(mlonp,precdatp,mloni)
    E0iprev=E0inext
    E0inext=reshape(parami,[lx2,lx3])
  else     !source data is fully 2D
    parami=interp2(mlonp,mlatp,Qp,mloni,mlati)     !interp to temp var.
    Qiprev=Qinext                       !save new pervious
    Qinext=reshape(parami,[lx2,lx3])    !overwrite next with new interpolated input

    parami=interp2(mlonp,mlatp,E0p,mloni,mlati)
    E0iprev=E0inext
    E0inext=reshape(parami,[lx2,lx3])
  end if
  if (myid==lid/2) then
    if (debug) print *, 'Min/max values for Qi:  ',minval(Qinext),maxval(Qinext)
    if (debug) print *, 'Min/max values for E0i:  ',minval(E0inext),maxval(E0inext)
  end if


  !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
  tprev=tnext
  UTsecprev=UTsecnext
  ymdprev=ymdnext

  tnext=tprev+cfg%dtprec
  UTsecnext=UTsectmp
  ymdnext=ymdtmp
end if


!INTERPOLATE IN TIME (LINEAR)
do ix3=1,lx3
  do ix2=1,lx2
    slope(ix2,ix3)=(Qinext(ix2,ix3)-Qiprev(ix2,ix3))/(tnext-tprev)
    Qinow(ix2,ix3)=Qiprev(ix2,ix3)+slope(ix2,ix3)*(t+dt/2._wp-tprev)

    slope(ix2,ix3)=(E0inext(ix2,ix3)-E0iprev(ix2,ix3))/(tnext-tprev)
    E0inow(ix2,ix3)=E0iprev(ix2,ix3)+slope(ix2,ix3)*(t+dt/2._wp-tprev)
  end do
end do


!SOME BASIC DIAGNOSTICS
if (myid==lid/2 .and. debug) then
  print *, 'tprev,t,tnext:  ',tprev,t+dt/2._wp,tnext
  print *, 'Min/max values for Qinow:  ',minval(Qinow),maxval(Qinow)
  print *, 'Min/max values for E0inow:  ',minval(E0inow),maxval(E0inow)
  print *, 'Min/max values for Qiprev:  ',minval(Qiprev),maxval(Qiprev)
  print *, 'Min/max values for E0prev:  ',minval(E0iprev),maxval(E0iprev)
  print *, 'Min/max values for Qinext:  ',minval(Qinext),maxval(Qinext)
  print *, 'Min/max values for E0next:  ',minval(E0inext),maxval(E0inext)
end if


!ASSIGN RESULTS OF INTERPOLATION TO OUTPUT VARIABLES
!background precipitation
W0pk=cfg%W0BG
PhiWpk=cfg%PhiWBG
do ix3=1,lx3
  do ix2=1,lx2
    W0(ix2,ix3,1)=W0pk
    PhiWmWm2(ix2,ix3,1)=PhiWpk
  end do
end do

!disturbance precipitation
W0(:,:,2)=E0inow
PhiWmWm2(:,:,2)=Qinow
end subroutine precipBCs_fileinput


subroutine clear_precip_fileinput()

if(allocated(mlonp)) then
  deallocate(mlonp,mlatp,Qp,E0p,mloni,mlati,Qiprev,Qinext,E0iprev,E0inext)
  if (allocated(precdatp)) then
    deallocate(precdatp)
  end if
end if

end subroutine clear_precip_fileinput


subroutine precipBCs(t,x,cfg,W0,PhiWmWm2)

!------------------------------------------------------------
!-------LOAD UP ARRAYS CONTAINING TOP BOUNDARY CHAR. ENERGY
!-------AND TOTAL ENERGY FLUX.  GRID VARIABLES INCLUDE
!-------GHOST CELLS
!------------------------------------------------------------

real(wp), intent(in) :: t
type(curvmesh), intent(in) :: x
type(gemini_cfg), intent(in) :: cfg
real(wp), dimension(:,:,:), intent(out) :: W0,PhiWmWm2

real(wp) :: W0pk,PhiWpk,meanW0x3,meanPhiWx3,sigW0x3,sigPhiWx3
real(wp) :: sigx2,meanx3,sigx3,x30amp,varc,meanx2,x2enve,sigt,meant
integer :: ix2,ix3,iprec,lx2,lx3,lprec


lx2=size(W0,1)
lx3=size(W0,2)
lprec=size(W0,3)    !assumed to be 2 in this subroutine


!BACKGROUND PRECIPITATION
W0pk = cfg%W0BG
PhiWpk=cfg%PhiWBG
!PhiWpk = 1e-3_wp
do ix3=1,lx3
  do ix2=1,lx2
    W0(ix2,ix3,1)=W0pk
    PhiWmWm2(ix2,ix3,1)=PhiWpk
  end do
end do


!PARAMETERS FOR DISTURBANCE PRECIPITATION
W0pk = 100
!      sigW0x3=100e3_wp
!      meanW0x3=0
PhiWpk = 0
!      PhiWpk=1e-5_wp    !successful grad-drift attempts
!      PhiWpk=1e-4_wp   !Swoboda blur testing
!      PhiWpk=0.05_wp    !testing of convergent Hall drifts
!      PhiWpk=5._wp
!      sigPhiWx3=100e3_wp
!      meanPhiWx3=0

!      W0pk=0.3e3_wp
!      sigW0x3=100e3_wp
!      meanW0x3=0
!      PhiWpk=2._wp
!      sigPhiWx3=100e3_wp
!      meanPhiWx3=0

sigx2 = 50e3_wp
meanx2 = 0
!    sigx3=10e3_wp
sigx3 = 25e3_wp
meant = 900
sigt = 450
x30amp= 0
varc = 200

!DISTURBANCE ELECTRON PRECIPITATION PATTERN
do ix3=1,lx3
  do ix2=1,lx2
    W0(ix2,ix3,2) = W0pk
    PhiWmWm2(ix2,ix3,2) = PhiWpk
  end do
end do

end subroutine precipBCs


end module precipBCs_mod
