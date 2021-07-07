module potentialBCs_mumps

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use mpimod, only : mpi_cfg
use phys_consts, only: wp, pi, Re, debug
use grid, only: lx1, lx2, lx2all, lx3all, gridflag
use meshobj, only: curvmesh
use interpolation, only : interp1,interp2
use timeutils, only : dateinc, date_filename, find_lastdate
use reader, only : get_grid2, get_simsize2, get_Efield
use config, only: gemini_cfg

implicit none (type, external)
private
public :: potentialbcs2D, potentialbcs2D_fileinput, clear_potential_fileinput, init_Efieldinput, &
            compute_rootBGEfields

!ALL OF THE FOLLOWING MODULE-SCOPE ARRAYS ARE USED FOR INTERPOLATING electric field INPUT FILES (IF USED)
!It should be noted that all of these will eventually be fullgrid variables since only root does this...
real(wp), dimension(:), allocatable :: mlonp
real(wp), dimension(:), allocatable :: mlatp    !coordinates of electric field data
integer :: llon,llat

real(wp), dimension(:,:), allocatable :: E0xp,E0yp    !x (lon.) and y (lat.) components of the electric field
real(wp), dimension(:,:), allocatable :: Vminx1p,Vmaxx1p
real(wp), dimension(:), allocatable :: Vminx2pslice,Vmaxx2pslice    !only slices because field lines (x1-dimension) should be equipotentials
real(wp), dimension(:), allocatable :: Vminx3pslice,Vmaxx3pslice
real(wp), dimension(:), allocatable :: Edatp    !needed when a 1D interpolation is to be done, i.e. when there is 1D sourde data

real(wp), dimension(:), allocatable :: mloni    !flat list of mlat,mlon locations on grid that we need to interpolate onto
real(wp), dimension(:), allocatable :: mlati

real(wp), dimension(:,:), allocatable :: E0xiprev,E0xinext,E0yiprev,E0yinext    !fields interpolated spatially
real(wp), dimension(:,:), allocatable :: Vminx1iprev,Vminx1inext,Vmaxx1iprev,Vmaxx1inext
real(wp), dimension(:), allocatable :: Vminx2isprev,Vminx2isnext,Vmaxx2isprev,Vmaxx2isnext
real(wp), dimension(:), allocatable :: Vminx3isprev,Vminx3isnext,Vmaxx3isprev,Vmaxx3isnext

real(wp), dimension(:,:), allocatable :: E0xinow,E0yinow

integer, dimension(3) :: ymdprev,ymdnext   !dates for interpolated data
real(wp) :: UTsecprev,UTsecnext
real(wp) :: tprev=0._wp,tnext=0._wp

integer :: ix1ref,ix2ref,ix3ref
!! reference location along field line closest to reference point of input data (300 km alt. at the grid center)

integer, private :: flagdirich_state
!! NOTE: holds state of flagdirich between calls, do not delete!

integer, private :: ix1eq=-1
!! index for the equatorial location in terms of index into the x%x1 array...

contains


subroutine init_Efieldinput(dt,t,cfg,ymd,UTsec,x)
  !> Initialize variables to hold electric field input file data, can be called by any worker but only root does anything

  !! FIXME:  This routine was extremely confusing to get right; much of that could probably be alleviated by not having t and ymd,UTsec
  !   assigned independently...

  real(wp), intent(in) :: dt,t
  type(gemini_cfg), intent(in) :: cfg
  integer, dimension(3), intent(in) :: ymd
  real(wp), intent(in) :: UTsec
  class(curvmesh), intent(in) :: x

  !! only root needs to allocate these
  real(wp), dimension(:,:), allocatable :: Vminx1,Vmaxx1
  real(wp), dimension(:,:), allocatable :: Vminx2,Vmaxx2
  real(wp), dimension(:,:), allocatable :: Vminx3,Vmaxx3
  real(wp), dimension(:,:,:), allocatable :: E01all,E02all,E03all
  integer :: flagdirich
  integer, dimension(3) :: ymdtmp
  real(wp) :: UTsectmp


  !> initializes the auroral electric field/current and particle inputs to read in a file corresponding to the first time step
  if (mpi_cfg%myid==0 .and. cfg%flagE0file==1) then    !only root needs these...
    allocate(Vminx1(1:x%lx2all,1:x%lx3all),Vmaxx1(1:x%lx2all,1:x%lx3all))
    allocate(Vminx2(1:x%lx1,1:x%lx3all),Vmaxx2(1:x%lx1,1:x%lx3all))
    allocate(Vminx3(1:x%lx1,1:x%lx2all),Vmaxx3(1:x%lx1,1:x%lx2all))
    allocate(E01all(1:x%lx1,1:x%lx2all,1:x%lx3all),E02all(1:x%lx1,1:x%lx2all,1:x%lx3all), &
                E03all(1:x%lx1,1:x%lx2all,1:x%lx3all))    !background fields

    !! find the last input data preceding the milestone/initial condition that we start with
    call find_lastdate(cfg%ymd0,cfg%UTsec0,ymd,UTsec,cfg%dtE0,ymdtmp,UTsectmp)

    print '(A,I4,I2.2,I2.2,F10.3)', 'init_Efieldinput: Prime electric field input files: ymd,utsec: ',ymdtmp,UTsectmp
    !! back up by one dt for each input so that we get a next file that corresponds to the first frame before
    !   this simulation beginning
    ! for first call the dummary argument must contain the input file time two levels back
    tprev=UTsectmp-UTsec-2._wp*cfg%dtE0
    tnext=tprev+cfg%dtE0
    call potentialBCs2D_fileinput(dt,tnext+dt/2._wp,ymdtmp,UTsectmp-cfg%dtE0,cfg,x,Vminx1, &
                                      Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                      Vmaxx3,E01all,E02all,E03all,flagdirich)    ! t input only needs to be less than zero...

    !! now load first, next frame of input corresponding to the initial time step.  Time input just needs to be negative to
    !   trigger a file read.  This might correspond to 0-dt/2???
    print*, 'init_Efieldinput: load next file for electric field input',tnext
    call potentialBCs2D_fileinput(dt,0._wp,ymdtmp,UTsectmp,cfg,x, &
                                      Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                      Vmaxx3,E01all,E02all,E03all,flagdirich)

    deallocate(Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3,E01all,E02all,E03all)
  end if
end subroutine init_Efieldinput


subroutine potentialBCs2D_fileinput(dt,t,ymd,UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                  Vmaxx3,E01all,E02all,E03all,flagdirich)
  !! A FILE INPUT BASED BOUNDARY CONDITIONS FOR ELECTRIC POTENTIAL OR
  !! FIELD-ALIGNED CURRENT.
  !! NOTE: THIS IS ONLY CALLED BY THE ROOT PROCESS

  real(wp), intent(in) :: dt
  real(wp), intent(in) :: t
  integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
  real(wp), intent(in) :: UTsec

  type(gemini_cfg), intent(in) :: cfg
  class(curvmesh), intent(in) :: x

  real(wp), dimension(:,:), intent(inout), target :: Vminx1,Vmaxx1
  !! intent(out)
  real(wp), dimension(:,:), intent(inout) :: Vminx2,Vmaxx2
  !! intent(out)
  real(wp), dimension(:,:), intent(inout) :: Vminx3,Vmaxx3
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: E01all,E02all,E03all
  !! intent(out)
  integer, intent(out) :: flagdirich

  real(wp) :: UTsectmp
  integer, dimension(3) :: ymdtmp

  real(wp), dimension(lx2all*lx3all) :: parami
  real(wp), dimension(lx2all,lx3all) :: parami2D
  real(wp), dimension(lx2all) :: parami2    !interpolated parameter with size of lx2
  real(wp), dimension(lx3all) :: parami3
  real(wp), dimension(lx2all,lx3all) :: Vminx1inow,Vmaxx1inow
  real(wp), dimension(lx3all) :: Vminx2isnow,Vmaxx2isnow
  real(wp), dimension(lx2all) :: Vminx3isnow,Vmaxx3isnow
  real(wp) :: slope

  integer :: ix1,ix2,ix3,iid,iflat,ios    !grid sizes are borrowed from grid module


  !> COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
  E01all = 0
  !! do not allow a background parallel field


  !! Note that the t<0 situation is going to cause the code to assume that ymd,UTsec have exactly
  !   the date/time of the input file that will be loaded...


  !FILE INPUT FOR THE PERPENDICULAR COMPONENTS OF THE ELECTRIC FIELD (ZONAL - X2, MERIDIONAL - X3)
  if(t + dt / 2._wp >= tnext .or. t < 0) then    !need to load a new file
    if ( .not. allocated(mlonp)) then    !need to read in the grid data from input file
      !! This is awful but it allows the initialization to get the correct files loaded for a restart...
      ymdprev=ymd
      UTsecprev=UTsec
      ymdnext=ymdprev
      UTsecnext=UTsecprev

      call get_simsize2(cfg%E0dir, llon=llon, llat=llat)

      if (debug) print '(A,2I6)', 'Electric field data has llon,llat size:  ',llon,llat
      if (llon < 1 .or. llat < 1) error stop 'potentialBCs_mumps: grid size must be strictly positive'
      allocate(mlonp(llon),mlatp(llat))
      !! bit of code duplication with worker code block below...


      !> IF WE HAVE SINGLETON DIMENSION, ALLOCATE SOME SPACE FOR TEMP ARRAY INPUT TO INTERP1
      if (llon==1) then
        allocate(Edatp(llat))
      elseif (llat==1) then
        allocate(Edatp(llon))
      end if


      !> NOW READ THE GRID
      call get_grid2(cfg%E0dir, mlonp, mlatp)

      if (debug) print *, 'Electric field data has mlon,mlat extent:', &
                minval(mlonp(:)), maxval(mlonp(:)), minval(mlatp(:)), maxval(mlatp(:))
      if(.not. all(ieee_is_finite(mlonp))) error stop 'mlon must be finite'
      if(.not. all(ieee_is_finite(mlatp))) error stop 'mlat must be finite'
      !> SPACE TO STORE INPUT DATA
      allocate(E0xp(llon,llat),E0yp(llon,llat))
      allocate(Vminx1p(llon,llat),Vmaxx1p(llon,llat))
      allocate(Vminx2pslice(llat),Vmaxx2pslice(llat))
      allocate(Vminx3pslice(llon),Vmaxx3pslice(llon))
      allocate(E0xiprev(lx2all,lx3all),E0xinext(lx2all,lx3all),E0yiprev(lx2all,lx3all),E0yinext(lx2all,lx3all))
      allocate(E0xinow(lx2all,lx3all),E0yinow(lx2all,lx3all))
      allocate(Vminx1iprev(lx2all,lx3all),Vminx1inext(lx2all,lx3all), &
               Vmaxx1iprev(lx2all,lx3all),Vmaxx1inext(lx2all,lx3all), &
               Vminx2isprev(lx3all),Vminx2isnext(lx3all),Vmaxx2isprev(lx3all),Vmaxx2isnext(lx3all), &
               Vminx3isprev(lx2all),Vminx3isnext(lx2all),Vmaxx3isprev(lx2all),Vmaxx3isnext(lx2all))

      E0xiprev = 0; E0yiprev = 0; E0xinext = 0; E0yinext = 0;
      E0xinow=0; E0yinow=0;
      !! these need to be initialized so that something sensible happens at the beginning
      Vminx1iprev = 0; Vmaxx1iprev = 0; Vminx1inext = 0; Vmaxx1inext = 0;
      Vminx2isprev = 0; Vmaxx2isprev = 0; Vminx2isnext = 0; Vmaxx2isnext = 0;
      Vminx3isprev = 0; Vmaxx3isprev = 0; Vminx3isnext = 0; Vmaxx3isnext = 0;


      !ALL PROCESSES NEED TO DEFINE THE POINTS THAT THEY WILL BE INTERPOLATING ONTO
      if (lx2all > 1 .and. lx3all>1) then ! 3D sim
        ix2ref = lx2all/2      !note integer division
        ix3ref = lx3all/2
      else if (lx2all==1 .and. lx3all>1) then
        ix2ref = 1
        ix3ref=lx3all/2
      else if (lx2all>1 .and. lx3all==1) then
        ix2ref=lx2all/2
        ix3ref=1
      else
        error stop 'Unable to orient boundary conditions for electric potential'
      endif


      !! by default the code uses 300km altitude as a reference location, using the center x2,x3 point
      ix1ref=minloc(abs(x%rall(:,ix2ref,ix3ref)-Re-300e3_wp),1)
      allocate(mloni(lx2all*lx3all),mlati(lx2all*lx3all))
      do ix3=1,lx3all
        do ix2=1,lx2all
          iflat=(ix3-1)*lx2all+ix2
          mlati(iflat)=90-x%thetaall(ix1ref,ix2,ix3)*180/pi
          mloni(iflat)=x%phiall(ix1ref,ix2,ix3)*180/pi
        end do
      end do
      if (debug) print '(A,4F7.2)', 'Grid has mlon,mlat range:  ',minval(mloni),maxval(mloni),minval(mlati),maxval(mlati)
      if (debug) print *, 'Grid has size:  ',iflat
    end if


    !> GRID INFORMATION EXISTS AT THIS POINT SO START READING IN PRECIP DATA
    !> read in the data from file
    if (debug) print *,'potentialBCs_mumps.f90:potentialBCs2D_fileinput: tprev,tnow,tnext:  ',tprev,t+dt/2,tnext
    ymdtmp=ymdnext
    UTsectmp=UTsecnext
    call dateinc(cfg%dtE0,ymdtmp,UTsectmp)
    !! get the date for "next" params

    call get_Efield(date_filename(cfg%E0dir, ymdtmp, UTsectmp), &
      flagdirich_state,E0xp,E0yp,Vminx1p,Vmaxx1p,&
      Vminx2pslice,Vmaxx2pslice,Vminx3pslice,Vmaxx3pslice)

    if (debug) then
      print *, 'Min/max values for E0xp:  ',minval(E0xp),maxval(E0xp)
      print *, 'Min/max values for E0yp:  ',minval(E0yp),maxval(E0yp)
      print *, 'Min/max values for Vminx1p:  ',minval(Vminx1p),maxval(Vminx1p)
      print *, 'Min/max values for Vmaxx1p:  ',minval(Vmaxx1p),maxval(Vmaxx1p)
      print *, 'Min/max values for Vminx2pslice:  ',minval(Vminx2pslice),maxval(Vminx2pslice)
      print *, 'Min/max values for Vmaxx2pslice:  ',minval(Vmaxx2pslice),maxval(Vmaxx2pslice)
      print *, 'Min/max values for Vminx3pslice:  ',minval(Vminx3pslice),maxval(Vminx3pslice)
      print *, 'Min/max values for Vmaxx3pslice:  ',minval(Vmaxx3pslice),maxval(Vmaxx3pslice)
    endif

    if (.not. all(ieee_is_finite(E0xp))) error stop 'E0xp: non-finite value(s)'
    if (.not. all(ieee_is_finite(E0yp))) error stop 'E0yp: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vminx1p))) error stop 'Vminx1p: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vmaxx1p))) error stop 'Vmaxx1p: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vminx2pslice))) error stop 'Vminx2pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vmaxx2pslice))) error stop 'Vmaxx2pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vminx3pslice))) error stop 'Vminx3pslice: non-finite value(s)'
    if (.not. all(ieee_is_finite(Vmaxx3pslice))) error stop 'Vmaxx3pslice: non-finite value(s)'

    !ALL WORKERS DO SPATIAL INTERPOLATION TO THEIR SPECIFIC GRID SITES
    if (debug) print *, 'Initiating electric field boundary condition spatial interpolations for date:  ',ymdtmp,' ',UTsectmp
    if (llon==1) then
      ! source data has singleton dimension in longitude
      if (debug) print *, 'Singleton longitude dimension detected; interpolating in latitude...'
      Edatp=E0xp(1,:)
      parami=interp1(mlatp,Edatp,mlati)
      !! will work even for 2D grids, just repeats the data in the lon direction
      E0xiprev=E0xinext
      E0xinext=reshape(parami,[lx2all,lx3all])

      Edatp=E0yp(1,:)
      parami=interp1(mlatp,Edatp,mlati)
      E0yiprev=E0yinext
      E0yinext=reshape(parami,[lx2all,lx3all])

      Edatp=Vminx1p(1,:)
      !! both min and max need to be read in from file and interpolated
      parami=interp1(mlatp,Edatp,mlati)
      Vminx1iprev=Vminx1inext
      Vminx1inext=reshape(parami,[lx2all,lx3all])

      Edatp=Vmaxx1p(1,:)
      parami=interp1(mlatp,Edatp,mlati)
      Vmaxx1iprev=Vmaxx1inext
      Vmaxx1inext=reshape(parami,[lx2all,lx3all])

      !! Note: for 2D simulations we don't use Vmaxx2p, etc. data read in from the input file - these BC's will be set later
    elseif (llat==1) then
      if (debug) print *, 'Singleton latitude dimension detected; interpolating in longitude...'
      Edatp=E0xp(:,1)
      parami=interp1(mlonp,Edatp,mloni)
      E0xiprev=E0xinext
      E0xinext=reshape(parami,[lx2all,lx3all])

      Edatp=E0yp(:,1)
      parami=interp1(mlonp,Edatp,mloni)
      E0yiprev=E0yinext
      E0yinext=reshape(parami,[lx2all,lx3all])

      Edatp=Vminx1p(:,1)
      parami=interp1(mlonp,Edatp,mloni)
      Vminx1iprev=Vminx1inext
      Vminx1inext=reshape(parami,[lx2all,lx3all])

      Edatp=Vmaxx1p(:,1)
      parami=interp1(mlonp,Edatp,mloni)
      Vmaxx1iprev=Vmaxx1inext
      Vmaxx1inext=reshape(parami,[lx2all,lx3all])
    else
      !! source data is 2D
      if (debug) print *, 'Executing full lat/lon interpolation...'
      parami=interp2(mlonp,mlatp,E0xp,mloni,mlati)     !interp to temp var.
      E0xiprev=E0xinext                       !save new previous
      E0xinext=reshape(parami,[lx2all,lx3all])    !overwrite next with new interpolated input

      parami=interp2(mlonp,mlatp,E0yp,mloni,mlati)
      E0yiprev=E0yinext
      E0yinext=reshape(parami,[lx2all,lx3all])

      parami=interp2(mlonp,mlatp,Vminx1p,mloni,mlati)
      Vminx1iprev=Vminx1inext
      Vminx1inext=reshape(parami,[lx2all,lx3all])

      parami=interp2(mlonp,mlatp,Vmaxx1p,mloni,mlati)
      Vmaxx1iprev=Vmaxx1inext
      Vmaxx1inext=reshape(parami,[lx2all,lx3all])

      !We need to interpolate the lateral boundaries in the direction of mlat
      parami=interp1(mlatp,Vminx2pslice,mlati)    !note mlati is a flat list of grid point lats, so need to reshape it
      Vminx2isprev=Vminx2isnext
      parami2D=reshape(parami,[lx2all,lx3all])
      parami3=parami2D(1,:)      !data should be constant across mlon, i.e. we're hoping the grid is plaid in mlat and mlon, otherwise not sure what to do here
      Vminx2isnext=parami3

      parami=interp1(mlatp,Vmaxx2pslice,mlati)
      Vmaxx2isprev=Vmaxx2isnext
      parami2D=reshape(parami,[lx2all,lx3all])
      parami3=parami2D(1,:)      !data should be constant across mlon...
      Vmaxx2isnext=parami3

      !now lateral interpolation in mlon
      parami=interp1(mlonp,Vminx3pslice,mloni)
      Vminx3isprev=Vminx3isnext
      parami2D=reshape(parami,[lx2all,lx3all])
      parami2=parami2D(:,1)
      Vminx3isnext=parami2

      parami=interp1(mlonp,Vmaxx3pslice,mloni)
      Vmaxx3isprev=Vmaxx3isnext
      parami2D=reshape(parami,[lx2all,lx3all])
      parami2=parami2D(:,1)
      Vmaxx3isnext=parami2
    end if

    if (debug) then
      print *, 'Min/max values for E0xi:  ',minval(E0xinext),maxval(E0xinext)
      print *, 'Min/max values for E0yi:  ',minval(E0yinext),maxval(E0yinext)
      print *, 'Min/max values for Vminx1i:  ',minval(Vminx1inext),maxval(Vminx1inext)
      print *, 'Min/max values for Vmaxx1i:  ',minval(Vmaxx1inext),maxval(Vmaxx1inext)

      if (llon/=1 .and. llat/=1) then
        print *, 'Min/max values for Vminx2i:  ',minval(Vminx2isnext),maxval(Vminx2isnext)
        print *, 'Min/max values for Vmaxx2i:  ',minval(Vmaxx2isnext),maxval(Vmaxx2isnext)
        print *, 'Min/max values for Vminx3i:  ',minval(Vminx3isnext),maxval(Vminx3isnext)
        print *, 'Min/max values for Vmaxx3i:  ',minval(Vmaxx3isnext),maxval(Vmaxx3isnext)
      end if
    end if


    !> UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
    tprev=tnext
    UTsecprev=UTsecnext
    ymdprev=ymdnext

    tnext=tprev+cfg%dtE0
    UTsecnext=UTsectmp
    ymdnext=ymdtmp
  end if


  !> INTERPOLATE IN TIME (LINEAR)
  flagdirich = flagdirich_state
  !! make sure to set solve type every time step, as it does not persiste between function calls
  if (debug) print *, 'Solve type: ',flagdirich
  do ix3=1,lx3all
    do ix2=1,lx2all
      slope=(E0xinext(ix2,ix3)-E0xiprev(ix2,ix3))/(tnext-tprev)
      E0xinow(ix2,ix3)=E0xiprev(ix2,ix3)+slope*(t+dt/2-tprev)

      slope=(E0yinext(ix2,ix3)-E0yiprev(ix2,ix3))/(tnext-tprev)
      E0yinow(ix2,ix3)=E0yiprev(ix2,ix3)+slope*(t+dt/2-tprev)

      slope=(Vminx1inext(ix2,ix3)-Vminx1iprev(ix2,ix3))/(tnext-tprev)
      Vminx1inow(ix2,ix3)=Vminx1iprev(ix2,ix3)+slope*(t+dt/2-tprev)

      slope=(Vmaxx1inext(ix2,ix3)-Vmaxx1iprev(ix2,ix3))/(tnext-tprev)
      Vmaxx1inow(ix2,ix3)=Vmaxx1iprev(ix2,ix3)+slope*(t+dt/2-tprev)
    end do
  end do
  if (lx2all/=1 .and. lx3all/=1) then
    !! full 3D grid need to also handle lateral boundaries; for 2D grids these will be set automatically later (or not used)
    do ix3=1,lx3all
      slope=(Vminx2isnext(ix3)-Vminx2isprev(ix3))/(tnext-tprev)
      Vminx2isnow(ix3)=Vminx2isprev(ix3)+slope*(t+dt/2-tprev)

      slope=(Vmaxx2isnext(ix3)-Vmaxx2isprev(ix3))/(tnext-tprev)
      Vmaxx2isnow(ix3)=Vmaxx2isprev(ix3)+slope*(t+dt/2-tprev)
    end do
    do ix2=1,lx2all
      slope=(Vminx3isnext(ix2)-Vminx3isprev(ix2))/(tnext-tprev)
      Vminx3isnow(ix2)=Vminx3isprev(ix2)+slope*(t+dt/2-tprev)

      slope=(Vmaxx3isnext(ix2)-Vmaxx3isprev(ix2))/(tnext-tprev)
      Vmaxx3isnow(ix2)=Vmaxx3isprev(ix2)+slope*(t+dt/2-tprev)
    end do
  end if


  !> SOME BASIC DIAGNOSTICS
  if(debug) then
    print *, 'tprev,t,tnext:  ',tprev,t+dt/2,tnext
    print *, 'Min/max values for E0xinow:  ',minval(E0xinow),maxval(E0xinow)
    print *, 'Min/max values for E0yinow:  ',minval(E0yinow),maxval(E0yinow)
    print *, 'Min/max values for Vminx1inow:  ',minval(Vminx1inow),maxval(Vminx1inow)
    print *, 'Min/max values for Vmaxx1inow:  ',minval(Vmaxx1inow),maxval(Vmaxx1inow)

    if (llon/=1 .and. llat/=1) then
      print *, 'Min/max values for Vminx2inow:  ',minval(Vminx2isnow),maxval(Vminx2isnow)
      print *, 'Min/max values for Vmaxx2inow:  ',minval(Vmaxx2isnow),maxval(Vmaxx2isnow)
      print *, 'Min/max values for Vminx3inow:  ',minval(Vminx3isnow),maxval(Vminx3isnow)
      print *, 'Min/max values for Vmaxx3inow:  ',minval(Vmaxx3isnow),maxval(Vmaxx3isnow)
    end if
  endif


  call compute_rootBGEfields(x,E02all,E03all)


  !> NOW THE BOUNDARY CONDITIONS
  do ix3=1,lx3all
    do ix2=1,lx2all
      Vminx1(ix2,ix3)=Vminx1inow(ix2,ix3)
      Vmaxx1(ix2,ix3)=Vmaxx1inow(ix2,ix3)
    end do
  end do


  ! This forces certain patterns to the boundary conditions to make sure solvers don't get garbage data...
  if (lx2all/=1 .and. lx3all/=1) then
    !! full 3D grid
    do ix3=1,lx3all
      do ix1=1,lx1
        Vminx2(ix1,ix3)=Vminx2isnow(ix3)
        Vmaxx2(ix1,ix3)=Vmaxx2isnow(ix3)
      end do
    end do

    do ix2=1,lx2all
      do ix1=1,lx1
        Vminx3(ix1,ix2)=Vminx3isnow(ix2)
        Vmaxx3(ix1,ix2)=Vmaxx3isnow(ix2)
      end do
    end do
  else
    !! some type of 2D grid, lateral boundary will be overwritten
    Vminx2 = 0
    Vmaxx2 = 0
    if (flagdirich==1) then
      !! Dirichlet:  needs to be the same as the physical top corner grid points
      if (gridflag/=1) then   !non-inverted so logical end is the max alt.
        if (lx2all==1) then
          do ix1=1,lx1
            Vminx3(ix1,:)=Vmaxx1(:,1)
            Vmaxx3(ix1,:)=Vmaxx1(:,lx3all)
          end do
          Vminx2=0
          Vmaxx2=0
        else
          do ix1=1,lx1
            Vminx2(ix1,:)=Vmaxx1(1,:)
            Vmaxx2(ix1,:)=Vmaxx1(lx2all,:)
          end do
          Vminx3=0
          Vmaxx3=0
        end if
      else                    !inverted so logical beginning is max alt.
        if (lx2all==1) then
          do ix1=1,lx1
            Vminx3(ix1,:)=Vminx1(:,1)
            Vmaxx3(ix1,:)=Vminx1(:,lx3all)
          end do
          Vminx2=0
          Vmaxx2=0
        else
          do ix1=1,lx1
            Vminx2(ix1,:)=Vminx1(1,:)
            Vmaxx2(ix1,:)=Vminx1(lx2all,:)
          end do
          Vminx3=0
          Vmaxx3=0
        end if
      end if
    else
      !! Neumann in x1:  sides are grounded...
      if (lx2all==1) then
        Vminx3 = 0
        Vmaxx3 = 0
      else
        Vminx2=0
        Vmaxx2=0
      end if
    end if
  end if
end subroutine potentialBCs2D_fileinput


subroutine compute_rootBGEfields(x,E02all,E03all)
  !> Returns a background electric field calculation for use by external program units.
  !   This requires that all necessary files, etc. have already been loaded into module
  !   variables.  This is only to be called by a root process as it deals with fullgrid
  !   data.  An interface for workers and root is in the top-level potential module.  This
  !   particular bit of code is needed both when setting boundary conditions and also when
  !   initializing background electric field; hence it is a subroutine as opposed to block of code

  class(curvmesh), intent(in) :: x
  real(wp), dimension(:,:,:), intent(inout) :: E02all,E03all
  !! intent(out)

  integer :: ix1,ix2,ix3
  real(wp) :: h2ref,h3ref


  !! the only danger here is that this routine could be called before any module data are loaded
  !   so check just to make sure it isn't being misused in this way
  if (.not. allocated(E0xinow)) error stop  &
        'potentialBCs:compute_rootBGEfields is trying to access unallocated module data'

  !! scale electric fields at some reference point into the full grid
  do ix3=1,lx3all
    do ix2=1,lx2all
      h2ref=x%h2all(ix1ref,ix2,ix3)
      !! define a reference metric factor for a given field line
      h3ref=x%h3all(ix1ref,ix2,ix3)
      do ix1=1,lx1
        E02all(ix1,ix2,ix3)=E0xinow(ix2,ix3)*h2ref/x%h2all(ix1,ix2,ix3)
        E03all(ix1,ix2,ix3)=E0yinow(ix2,ix3)*h3ref/x%h3all(ix1,ix2,ix3)
      end do
    end do
  end do

  end subroutine compute_rootBGEfields


  subroutine clear_potential_fileinput()

  if(allocated(mlonp)) then
    deallocate(mlonp,mlatp,mloni,mlati,E0xp,E0yp,Vminx1p,Vmaxx1p)
    deallocate(Vminx2pslice,Vmaxx2pslice)
    deallocate(Vminx3pslice,Vmaxx3pslice)
    if (allocated(Edatp)) then
      deallocate(Edatp)
    end if
    deallocate(E0xiprev,E0xinext,E0yiprev,E0yinext)
    deallocate(Vminx1iprev,Vminx1inext,Vmaxx1iprev,Vmaxx1inext)
    deallocate(Vminx2isprev,Vminx2isnext,Vmaxx2isprev,Vmaxx2isnext)
    deallocate(Vminx3isprev,Vminx3isnext,Vmaxx3isprev,Vmaxx3isnext)
  end if
end subroutine clear_potential_fileinput


subroutine potentialBCs2D(UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                      Vmaxx3,E01all,E02all,E03all,flagdirich)
  ! This is a default routine for setting electromagnetic boundary conditions in cases where user file input is not specified.  It also computes the equatorial vertical drift for the EIA if requested by the user. This routine *could* be modified to hard-code specific conditions in if needed but we really recommend using file input for that.

  real(wp), intent(in) :: UTsec

  type(gemini_cfg), intent(in) :: cfg
  class(curvmesh), intent(in) :: x

  real(wp), dimension(:,:), intent(inout), target :: Vminx1,Vmaxx1
  !! intent(out)
  real(wp), dimension(:,:), intent(inout) :: Vminx2,Vmaxx2
  !! intent(out)
  real(wp), dimension(:,:), intent(inout) :: Vminx3,Vmaxx3
  !! intent(out)
  real(wp), dimension(:,:,:), intent(inout) :: E01all,E02all,E03all
  !! intent(out)
  integer, intent(out) :: flagdirich

  real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field

  real(wp) :: Phipk
  integer :: ix1,ix2,ix3    !grid sizes are borrowed from grid module
  integer :: im
  !    integer, parameter :: lmodes=8; this type of thing done from input scripts now...
  real(wp) :: phase
  real(wp), dimension(1:size(Vmaxx1,1)) :: x3dev
  real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field

  real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt

  real(wp) :: vamp,LThrs,veltime,z,glonmer


  !CALCULATE/SET TOP BOUNDARY CONDITIONS
  sigx2=1/20._wp*(x%x2all(lx2all)-x%x2all(1))
  meanx2=0.5_wp*(x%x2all(1)+x%x2all(lx2all))
  sigx3=1/20._wp*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
  meanx3=0.5_wp*(x%x3all(1)+x%x3all(lx3all))


  ! FIXME: the pointer swapping to deal with top vs. bottom here is confusing/superfluous; it may be better simply to have the input preparation scripts assign things accordingly, which is what we will do from now on.  For this routine right now it doesn't matter since both just zeroed out anyway...
  if (gridflag/=2) then
    Vtopalt=>Vminx1
    Vbotalt=>Vmaxx1
  else
    Vtopalt=>Vmaxx1
    Vbotalt=>Vminx1
  end if

  Phipk = 0      !pk current density
  flagdirich = 0    !Neumann conditions
  do ix3=1,lx3all
    do ix2=1,lx2all
      Vtopalt(ix2,ix3) = 0
    end do
  end do


  !SOME USER INFO
  if (debug) print *, 'At time (UT seconds):  ',UTsec,'  Max FAC set to be:  ',maxval(abs(Vtopalt))


  !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
  Vbotalt = 0   !since we need to have no current through bottom boundary
  Vminx2 = 0
  Vmaxx2 = 0
  Vminx3 = 0
  Vmaxx3 = 0


  !PI's EIA code COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
  if (cfg%flagEIA) then
    if (ix1eq<=0) then   !recompute the position of the equator in terms of the x1 variable
      ix1eq=minloc(abs(x%x1),1)    !equator location is that closest to zero in the x1 (q) variable
      if (debug) print*, 'equator ix1:  ',ix1eq,x%x1(ix1eq)
    end if

    vamp=cfg%v0equator    !amplitude of vertical drift at equator from input config.nml file

    E01all=0
    E02all=0
    do ix2=1,lx2all    !for a swapped grid this is longitude
      !for each meridional slice define a local time
      glonmer=x%glonall(ix1eq,ix2,lx3all/2)     !just use halfway up in altitude at the magnetic equator
      do while (glonmer<0)
        glonmer=glonmer+360
      end do

      LThrs=UTsec/3600+glonmer/15                 !Local time of center of meridian
      veltime = sin(2*pi*(LThrs-7)/24)    ! Huba's formulate for velocity amplitude vs. time

      do ix3=1,lx3all     !here this is L-shell
        z = x%altall(ix1eq,ix2,ix3)  !Current altitude of center of this flux tube
        do ix1=1,lx1
          if (z<=150e3_wp) then
            E03all(ix1,ix2,ix3) = 0
          elseif ((z>=150e3_wp) .and. (z<=300e3_wp)) then
            E03all(ix1,ix2,ix3) = -(veltime*vamp*(z-150e3_wp)/150e3_wp)*x%Bmagall(ix1eq,ix2,ix3)    !minus sign to deal with permuted dimensions
          elseif (z>300e3_wp) then
            E03all(ix1,ix2,ix3) = -veltime*vamp*x%Bmagall(ix1eq,ix2,ix3)
          end if
        end do
      end do
    end do

  !  print*, '  Applied EIA perturbation to background electric field...'
  !  print*, '    ',minval(E03all),maxval(E03all)
  else
    E01all = 0
    E02all = 0
    E03all = 0
  end if
end subroutine potentialBCs2D


impure elemental subroutine assert_file_exists(path)
  !! throw error if file does not exist
  !! this accommodates non-Fortran 2018 error stop with variable character

  character(*), intent(in) :: path
  logical :: exists

  inquire(file=path, exist=exists)

  if (.not.exists) then
    write(stderr,*) path // ' does not exist'
    error stop
  endif
end subroutine assert_file_exists

end module potentialBCs_mumps
