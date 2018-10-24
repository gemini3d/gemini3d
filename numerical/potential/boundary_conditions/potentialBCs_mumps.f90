module potentialBCs_mumps

use grid, only: curvmesh, lx1, lx2, lx3all, gridflag
use phys_consts
use calculus
use interpolation, only : interp1,interp2
use io, only : date_filename
use temporal, only : dateinc

implicit none

!ALL OF THE FOLLOWING MODULE-SCOPE ARRAYS ARE USED FOR INTERPOLATING PRECIPITATION INPUT FILES (IF USED)
real(wp), dimension(:), allocatable, private :: mlonp
real(wp), dimension(:), allocatable, private :: mlatp    !coordinates of electric field data
integer, private :: llon,llat

real(wp), dimension(:,:), allocatable, private :: E0xp,E0yp    !x (lon.) and y (lat.) components of the electric field
real(wp), dimension(:,:), allocatable, private :: Vminx1p,Vmaxx1p
real(wp), dimension(:), allocatable, private :: Vminx2pslice,Vmaxx2pslice    !only slices because field lines (x1-dimension) should be equipotentials
real(wp), dimension(:), allocatable, private :: Vminx3pslice,Vmaxx3pslice
real(wp), dimension(:), allocatable, private :: Edatp    !needed when a 1D interpolation is to be done, i.e. when there is 1D sourde data

real(wp), dimension(:), allocatable, private :: mloni    !flat list of mlat,mlon locations on grid that we need to interpolate onto
real(wp), dimension(:), allocatable, private :: mlati

real(wp), dimension(:,:), allocatable, private :: E0xiprev,E0xinext,E0yiprev,E0yinext    !fields interpolated spatially
real(wp), dimension(:,:), allocatable, private :: Vminx1iprev,Vminx1inext,Vmaxx1iprev,Vmaxx1inext
real(wp), dimension(:), allocatable, private :: Vminx2isprev,Vminx2isnext,Vmaxx2isprev,Vmaxx2isnext
real(wp), dimension(:), allocatable, private :: Vminx3isprev,Vminx3isnext,Vmaxx3isprev,Vmaxx3isnext

integer, dimension(3), private :: ymdprev,ymdnext   !dates for interpolated data
real(wp), private :: UTsecprev,UTsecnext
real(wp), private :: tprev,tnext

real(wp), private :: flagdirich_double

contains


  subroutine potentialBCs2D_fileinput(dt,dtE0,t,ymd,UTsec,E0dir,&
                                      x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                      Vmaxx3,E01all,E02all,E03all,flagdirich)

    !A FILE INPUT BASED BOUNDARY CONDITIONS FOR ELECTRIC POTENTIAL OR
    !FIELD-ALIGNED CURRENT.  NOTE THAT THIS IS ONLY CALLED BY THE ROOT
    !PROCESS!!!

    real(wp), intent(in) :: dt
    real(wp), intent(in) :: dtE0    !cadence at which we are reading in the E0 files
    real(wp), intent(in) :: t
    integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    character(*), intent(in) :: E0dir       !directory where data are kept

    type(curvmesh), intent(in) :: x

    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
    integer, intent(out) :: flagdirich

    character(:), allocatable :: fn1, fn2, fn3
    integer :: inunit

    real(wp) :: UTsectmp
    integer, dimension(3) :: ymdtmp

    real(wp), dimension(lx2*lx3all) :: parami
    real(wp), dimension(lx2,lx3all) :: parami2D
    real(wp), dimension(lx2) :: parami2    !interpolated parameter with size of lx2
    real(wp), dimension(lx3all) :: parami3
    real(wp), dimension(lx2,lx3all) :: E0xinow,E0yinow,Vminx1inow,Vmaxx1inow
    real(wp), dimension(lx3all) :: Vminx2isnow,Vmaxx2isnow
    real(wp), dimension(lx2) :: Vminx3isnow,Vmaxx3isnow
    real(wp) :: slope

    integer :: ix1,ix2,ix3,iid,iflat,ios    !grid sizes are borrowed from grid module


    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
    E01all=0d0    !do not allow a background parallel field


    !FILE INPUT FOR THE PERPENDICULAR COMPONENTS OF THE ELECTRIC FIELD (ZONAL - X2, MERIDIONAL - X3)
    if(t+dt/2d0>=tnext) then    !need to load a new file
      if ( .not. allocated(mlonp)) then    !need to read in the grid data from input file
        ymdprev=ymd
        UTsecprev=UTsec
        ymdnext=ymdprev
        UTsecnext=UTsecprev

        fn1 = trim(E0dir) // '/simsize.dat'
        print *, 'Inputting electric field data size from file:  ',fn1
        open(newunit=inunit,file=fn1,status='old',form='unformatted',access='stream')
        read(inunit) llon,llat
        close(inunit)
        print *, 'Electric field data has llon,llat size:  ',llon,llat
        allocate(mlonp(llon),mlatp(llat))    !bit of code duplication with worker code block below...


        !IF WE HAVE SINGLETON DIMENSION THEN ALLOCATE SOME SPACE FOR A TEMP
        !ARRAY FOR INPUTTING INTO INTERP1
        if (llon==1) then
          allocate(Edatp(llat))
        elseif (llat==1) then
          allocate(Edatp(llon))
        end if


        !NOW READ THE GRID
        fn2 = E0dir // '/simgrid.dat'
        print *, 'Inputting electric field grid from file:  ',fn2
        open(newunit=inunit,file=fn2,status='old',form='unformatted',access='stream')
        read(inunit) mlonp,mlatp
        close(inunit)
        print *, 'Electric field data has mlon,mlat extent:', &
                  minval(mlonp(:)), maxval(mlonp(:)), minval(mlatp(:)), maxval(mlatp(:))

        !SPACE TO STORE INPUT DATA
        allocate(E0xp(llon,llat),E0yp(llon,llat))
        allocate(Vminx1p(llon,llat),Vmaxx1p(llon,llat))
        allocate(Vminx2pslice(llat),Vmaxx2pslice(llat))
        allocate(Vminx3pslice(llon),Vmaxx3pslice(llon))
        allocate(E0xiprev(lx2,lx3all),E0xinext(lx2,lx3all),E0yiprev(lx2,lx3all),E0yinext(lx2,lx3all))
        allocate(Vminx1iprev(lx2,lx3all),Vminx1inext(lx2,lx3all), &
                 Vmaxx1iprev(lx2,lx3all),Vmaxx1inext(lx2,lx3all), &
                 Vminx2isprev(lx3all),Vminx2isnext(lx3all),Vmaxx2isprev(lx3all),Vmaxx2isnext(lx3all), &
                 Vminx3isprev(lx2),Vminx3isnext(lx2),Vmaxx3isprev(lx2),Vmaxx3isnext(lx2))

        E0xiprev=0d0; E0yiprev=0d0; E0xinext=0d0; E0yinext=0d0;     !these need to be initialized so that something sensible happens at the beginning
        Vminx1iprev=0d0; Vmaxx1iprev=0d0; Vminx1inext=0d0; Vmaxx1inext=0d0;
        Vminx2isprev=0d0; Vmaxx2isprev=0d0; Vminx2isnext=0d0; Vmaxx2isnext=0d0;
        Vminx3isprev=0d0; Vmaxx3isprev=0d0; Vminx3isnext=0d0; Vmaxx3isnext=0d0;


        !ALL PROCESSES NEED TO DEFINED THE OPINTS THAT THEY WILL BE INTERPOLATING ONTO
        allocate(mloni(lx2*lx3all),mlati(lx2*lx3all))
        do ix3=1,lx3all
          do ix2=1,lx2
            iflat=(ix3-1)*lx2+ix2
            mlati(iflat)=90d0-x%thetaall(lx1,ix2,ix3)*180d0/pi    !does rool even have this info full grid???
            mloni(iflat)=x%phiall(lx1,ix2,ix3)*180d0/pi
          end do
        end do
        write(*,*) 'Grid has mlon,mlat range:  ',minval(mloni),maxval(mloni),minval(mlati),maxval(mlati)
        write(*,*) 'Grid has size:  ',iflat
      end if


      !GRID INFORMATION EXISTS AT THIS POINT SO START READING IN PRECIP DATA
      !read in the data from file
      print *,'tprev,tnow,tnext:  ',tprev,t+dt/2d0,tnext
      ymdtmp=ymdnext
      UTsectmp=UTsecnext
      call dateinc(dtE0,ymdtmp,UTsectmp)    !get the date for "next" params
      fn3=date_filename(E0dir,ymdtmp,UTsectmp)     !form the standard data filename
      print *, 'Pulling electric field data from file:  ',fn3
      open(newunit=inunit, file=fn3, status='old',form='unformatted',access='stream',iostat=ios)
      if (ios==0) then    !successful read
        read(inunit) flagdirich_double
        read(inunit) E0xp,E0yp
        read(inunit) Vminx1p,Vmaxx1p    !background fields and top/bottom boundar conditions
        read(inunit) Vminx2pslice,Vmaxx2pslice    !these ohly used for 3D simulations
        read(inunit) Vminx3pslice,Vmaxx3pslice
      else      
        error stop 'Bad input file, cannot proceed'  ! per MZ Oct. 2018
        !just set everything to zero
        !write(*,*) 'Bad input file, setting everything to some default value...'
        !flagdirich=1    !to short-circuit solve...
        !E0xp=0d0; E0yp=0d0; Vminx1p=0d0; Vmaxx1p=0d0;
        !Vminx2pslice=0d0; Vmaxx2pslice=0d0; Vminx3pslice=0d0; Vmaxx3pslice=0d0;
      end if
      close(inunit)

      write(*,*) 'Min/max values for E0xp:  ',minval(pack(E0xp,.true.)),maxval(pack(E0xp,.true.))
      write(*,*) 'Min/max values for E0yp:  ',minval(pack(E0yp,.true.)),maxval(pack(E0yp,.true.))
      write(*,*) 'Min/max values for Vminx1p:  ',minval(pack(Vminx1p,.true.)),maxval(pack(Vminx1p,.true.))
      write(*,*) 'Min/max values for Vmaxx1p:  ',minval(pack(Vmaxx1p,.true.)),maxval(pack(Vmaxx1p,.true.))
      write(*,*) 'Min/max values for Vminx2pslice:  ',minval(pack(Vminx2pslice,.true.)),maxval(pack(Vminx2pslice,.true.))
      write(*,*) 'Min/max values for Vmaxx2pslice:  ',minval(pack(Vmaxx2pslice,.true.)),maxval(pack(Vmaxx2pslice,.true.))
      write(*,*) 'Min/max values for Vminx3pslice:  ',minval(pack(Vminx3pslice,.true.)),maxval(pack(Vminx3pslice,.true.))
      write(*,*) 'Min/max values for Vmaxx3pslice:  ',minval(pack(Vmaxx3pslice,.true.)),maxval(pack(Vmaxx3pslice,.true.))


      !ALL WORKERS DO SPATIAL INTERPOLATION TO THEIR SPECIFIC GRID SITES
      write(*,*) 'Initiating electric field boundary condition spatial interpolations for date:  ',ymdtmp,' ',UTsectmp
      if (llon==1) then    !source data has singleton dimension in longitude
        write(*,*) 'Singleton longitude dimension detected; interpolating in latitude...'
        Edatp=E0xp(1,:)
        parami=interp1(mlatp,Edatp,mlati)   !will work even for 2D grids, just repeats the data in the lon direction
        E0xiprev=E0xinext
        E0xinext=reshape(parami,[lx2,lx3all])

        Edatp=E0yp(1,:)
        parami=interp1(mlatp,Edatp,mlati)
        E0yiprev=E0yinext
        E0yinext=reshape(parami,[lx2,lx3all])

        Edatp=Vminx1p(1,:)          !both min and max need to be read in from file and interpolated
        parami=interp1(mlatp,Edatp,mlati)
        Vminx1iprev=Vminx1inext
        Vminx1inext=reshape(parami,[lx2,lx3all])

        Edatp=Vmaxx1p(1,:)
        parami=interp1(mlatp,Edatp,mlati)
        Vmaxx1iprev=Vmaxx1inext
        Vmaxx1inext=reshape(parami,[lx2,lx3all])

        !note that for 2D simulations we don't use Vmaxx2p, etc. data read in from the input file - these BC's will be set later
      elseif (llat==1) then
        write(*,*) 'Singleton latitude dimension detected; interpolating in longitude...'
        Edatp=E0xp(:,1)
        parami=interp1(mlonp,Edatp,mloni)
        E0xiprev=E0xinext
        E0xinext=reshape(parami,[lx2,lx3all])

        Edatp=E0yp(:,1)
        parami=interp1(mlonp,Edatp,mloni)
        E0yiprev=E0yinext
        E0yinext=reshape(parami,[lx2,lx3all])

        Edatp=Vminx1p(:,1)
        parami=interp1(mlonp,Edatp,mloni)
        Vminx1iprev=Vminx1inext
        Vminx1inext=reshape(parami,[lx2,lx3all])

        Edatp=Vmaxx1p(:,1)
        parami=interp1(mlonp,Edatp,mloni)
        Vmaxx1iprev=Vmaxx1inext
        Vmaxx1inext=reshape(parami,[lx2,lx3all])
      else    !source data is 2D
        write(*,*) 'Executing full lat/lon interpolation...'
        parami=interp2(mlonp,mlatp,E0xp,mloni,mlati)     !interp to temp var.
        E0xiprev=E0xinext                       !save new pervious
        E0xinext=reshape(parami,[lx2,lx3all])    !overwrite next with new interpolated input

        parami=interp2(mlonp,mlatp,E0yp,mloni,mlati)
        E0yiprev=E0yinext
        E0yinext=reshape(parami,[lx2,lx3all])

        parami=interp2(mlonp,mlatp,Vminx1p,mloni,mlati)
        Vminx1iprev=Vminx1inext
        Vminx1inext=reshape(parami,[lx2,lx3all])

        parami=interp2(mlonp,mlatp,Vmaxx1p,mloni,mlati)
        Vmaxx1iprev=Vmaxx1inext
        Vmaxx1inext=reshape(parami,[lx2,lx3all])

        !We need to interpolate the lateral boundaries in the direction of mlat
        parami=interp1(mlatp,Vminx2pslice,mlati)    !note mlati is a flat list of grid point lats, so need to reshape it
        Vminx2isprev=Vminx2isnext
        parami2D=reshape(parami,[lx2,lx3all])
        parami3=parami2D(1,:)      !data should be constant across mlon, i.e. we're hoping the grid is plaid in mlat and mlon, otherwise not sure what to do here
        Vminx2isnext=parami3

        parami=interp1(mlatp,Vmaxx2pslice,mlati)
        Vmaxx2isprev=Vmaxx2isnext
        parami2D=reshape(parami,[lx2,lx3all])
        parami3=parami2D(1,:)      !data should be constant across mlon...
        Vmaxx2isnext=parami3

        !now lateral interpolation in mlon
        parami=interp1(mlonp,Vminx3pslice,mloni)
        Vminx3isprev=Vminx3isnext
        parami2D=reshape(parami,[lx2,lx3all])
        parami2=parami2D(:,1)
        Vminx3isnext=parami2

        parami=interp1(mlonp,Vmaxx3pslice,mloni)
        Vmaxx3isprev=Vmaxx3isnext
        parami2D=reshape(parami,[lx2,lx3all])
        parami2=parami2D(:,1)
        Vmaxx3isnext=parami2
      end if

      write(*,*) 'Min/max values for E0xi:  ',minval(pack(E0xinext,.true.)),maxval(pack(E0xinext,.true.))
      write(*,*) 'Min/max values for E0yi:  ',minval(pack(E0yinext,.true.)),maxval(pack(E0yinext,.true.))
      write(*,*) 'Min/max values for Vminx1i:  ',minval(pack(Vminx1inext,.true.)),maxval(pack(Vminx1inext,.true.))
      write(*,*) 'Min/max values for Vmaxx1i:  ',minval(pack(Vmaxx1inext,.true.)),maxval(pack(Vmaxx1inext,.true.))
      if (llon/=1 .and. llat/=1) then
        write(*,*) 'Min/max values for Vminx2i:  ',minval(pack(Vminx2isnext,.true.)),maxval(pack(Vminx2isnext,.true.))
        write(*,*) 'Min/max values for Vmaxx2i:  ',minval(pack(Vmaxx2isnext,.true.)),maxval(pack(Vmaxx2isnext,.true.))
        write(*,*) 'Min/max values for Vminx3i:  ',minval(pack(Vminx3isnext,.true.)),maxval(pack(Vminx3isnext,.true.))
        write(*,*) 'Min/max values for Vmaxx3i:  ',minval(pack(Vmaxx3isnext,.true.)),maxval(pack(Vmaxx3isnext,.true.))
      end if


      !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
      tprev=tnext
      UTsecprev=UTsecnext
      ymdprev=ymdnext

      tnext=tprev+dtE0
      UTsecnext=UTsectmp
      ymdnext=ymdtmp
    end if


    !INTERPOLATE IN TIME (LINEAR)
    flagdirich=int(flagdirich_double,4)     !make sure to set solve type every time step, as it does not persiste between function calls
    write(*,*) 'Solve type: ',flagdirich 
    do ix3=1,lx3all
      do ix2=1,lx2
        slope=(E0xinext(ix2,ix3)-E0xiprev(ix2,ix3))/(tnext-tprev)
        E0xinow(ix2,ix3)=E0xiprev(ix2,ix3)+slope*(t+dt/2d0-tprev)

        slope=(E0yinext(ix2,ix3)-E0yiprev(ix2,ix3))/(tnext-tprev)
        E0yinow(ix2,ix3)=E0yiprev(ix2,ix3)+slope*(t+dt/2d0-tprev)

        slope=(Vminx1inext(ix2,ix3)-Vminx1iprev(ix2,ix3))/(tnext-tprev)
        Vminx1inow(ix2,ix3)=Vminx1iprev(ix2,ix3)+slope*(t+dt/2d0-tprev)

        slope=(Vmaxx1inext(ix2,ix3)-Vmaxx1iprev(ix2,ix3))/(tnext-tprev)
        Vmaxx1inow(ix2,ix3)=Vmaxx1iprev(ix2,ix3)+slope*(t+dt/2d0-tprev)
      end do
    end do
    if (lx2/=1 .and. lx3all/=1) then     !full 3D grid need to also handle lateral boundaries
      do ix3=1,lx3all
        slope=(Vminx2isnext(ix3)-Vminx2isprev(ix3))/(tnext-tprev)
        Vminx2isnow(ix3)=Vminx2isprev(ix3)+slope*(t+dt/2-tprev)

        slope=(Vmaxx2isnext(ix3)-Vmaxx2isprev(ix3))/(tnext-tprev)
        Vmaxx2isnow(ix3)=Vmaxx2isprev(ix3)+slope*(t+dt/2-tprev)
      end do
      do ix2=1,lx2
        slope=(Vminx3isnext(ix2)-Vminx3isprev(ix2))/(tnext-tprev)
        Vminx3isnow(ix2)=Vminx3isprev(ix2)+slope*(t+dt/2-tprev)

        slope=(Vmaxx3isnext(ix2)-Vmaxx3isprev(ix2))/(tnext-tprev)
        Vmaxx3isnow(ix2)=Vmaxx3isprev(ix2)+slope*(t+dt/2-tprev)
      end do
    end if


    !SOME BASIC DIAGNOSTICS
    write(*,*) 'tprev,t,tnext:  ',tprev,t+dt/2d0,tnext
    write(*,*) 'Min/max values for E0xinow:  ',minval(pack(E0xinow,.true.)),maxval(pack(E0xinow,.true.))
    write(*,*) 'Min/max values for E0yinow:  ',minval(pack(E0yinow,.true.)),maxval(pack(E0yinow,.true.))
    write(*,*) 'Min/max values for Vminx1inow:  ',minval(pack(Vminx1inow,.true.)),maxval(pack(Vminx1inow,.true.))
    write(*,*) 'Min/max values for Vmaxx1inow:  ',minval(pack(Vmaxx1inow,.true.)),maxval(pack(Vmaxx1inow,.true.))
    if (llon/=1 .and. llat/=1) then
      write(*,*) 'Min/max values for Vminx2inow:  ',minval(pack(Vminx2isnow,.true.)),maxval(pack(Vminx2isnow,.true.))
      write(*,*) 'Min/max values for Vmaxx2inow:  ',minval(pack(Vmaxx2isnow,.true.)),maxval(pack(Vmaxx2isnow,.true.))
      write(*,*) 'Min/max values for Vminx3inow:  ',minval(pack(Vminx3isnow,.true.)),maxval(pack(Vminx3isnow,.true.))
      write(*,*) 'Min/max values for Vmaxx3inow:  ',minval(pack(Vmaxx3isnow,.true.)),maxval(pack(Vmaxx3isnow,.true.))
    end if


    !LOAD POTENTIAL SOLVER INPUT ARRAYS
    do ix3=1,lx3all
      do ix2=1,lx2
        do ix1=1,lx1
          E02all(ix1,ix2,ix3)=E0xinow(ix2,ix3)
          E03all(ix1,ix2,ix3)=E0yinow(ix2,ix3)
        end do
      end do
    end do
    do ix3=1,lx3all
      do ix2=1,lx2
        Vminx1(ix2,ix3)=Vminx1inow(ix2,ix3)
        Vmaxx1(ix2,ix3)=Vmaxx1inow(ix2,ix3)
      end do
    end do


    !SET REMAINING BOUNDARY CONDITIONS BASED ON WHAT THE TOP IS.  IF WE HAVE A
    !3D GRID THE SIDES ARE GROUNDED AUTOMATICALLY, WHEREAS FOR 2D THEY ARE SET
    !TO TOP VALUE  IF DIRICHLET AND TO TOP VALUE IF DIRICHLET.
    if (lx2/=1 .and. lx3all/=1) then     !full 3D grid
!      Vminx2=0d0    !This actualy needs to be different for KHI
!      Vmaxx2=0d0
!      Vminx3=0d0
!      Vmaxx3=0d0
      do ix3=1,lx3all
        do ix1=1,lx1
          Vminx2(ix1,ix3)=Vminx2isnow(ix3)
          Vmaxx2(ix1,ix3)=Vmaxx2isnow(ix3)
        end do
      end do

      do ix2=1,lx2
        do ix1=1,lx1
          Vminx3(ix1,ix2)=Vminx3isnow(ix2)
          Vmaxx3(ix1,ix2)=Vmaxx3isnow(ix2)
        end do
      end do
    else    !some type of 2D grid, lateral boundary will be overwritten
      Vminx2=0d0
      Vmaxx2=0d0
      if (flagdirich==1) then    !Dirichlet:  needs to be the same as the top corner grid points
        do ix1=1,lx1
          Vminx3(ix1,:)=Vmaxx1(:,1)
          Vmaxx3(ix1,:)=Vmaxx1(:,lx3all)
        end do
      else    !Neumann in x1:  sides are grounded...
        Vminx3=0d0
        Vmaxx3=0d0
      end if
    end if

  end subroutine potentialBCs2D_fileinput


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


!  subroutine potentialBCs2D_fileinput(dt,dtE0,t,ymd,UTsec,E0dir,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !A FILE INPUT BASED BOUNDARY CONDITIONS FOR ELECTRIC POTENTIAL OR
!    !FIELD-ALIGNED CURRENT.  NOTE THAT THIS IS ONLY CALLED BY THE ROOT
!    !PROCESS!!!
!
!    real(wp), intent(in) :: dt
!    real(wp), intent(in) :: dtE0    !cadence at which we are reading in the E0 files
!    real(wp), intent(in) :: t
!    integer, dimension(3), intent(in) :: ymd    !date for which we wish to calculate perturbations
!    real(wp), intent(in) :: UTsec
!    character(*), intent(in) :: E0dir       !directory where data is kept
!
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    character(512) :: filename 
!    integer, parameter :: inunit=48
!
!    real(wp) :: UTsectmp
!    integer, dimension(3) :: ymdtmp
!
!    real(wp), dimension(lx2*lx3all) :: parami
!    real(wp), dimension(lx2,lx3all) :: slope,E0xinow,E0yinow    !these are for each worker
!    real(wp), dimension(:,:), allocatable :: E0xinowall,E0yinowall    !full grid copies for root to assemble
!
!    integer :: ix1,ix2,ix3,iid,iflat,ios    !grid sizes are borrowed from grid module
!
!
!
!    !FIRST BIT OF CODE RESETS THE BOUNDARIES AND PARALLEL COMPONENTS WHICH SHOULD BE FIXED AT 0D0
!    flagdirich=1    !Should be based on file input, as well.  This may also set some of the boundary conditions - refer to solver source code 
!
!
!    !EVERYONE IS GROUNDED TO ALLOW FOR BOTH BACKGROUND FIELDS AND FACS
!    Vminx1=0d0   !since we need to have no current through bottom boundary
!    Vmaxx1=0d0
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0    !do not allow a background parallel field
!
!
!    !FILE INPUT FOR THE PERPENDICULAR COMPONENTS OF THE ELECTRIC FIELD (ZONAL - X2, MERIDIONAL - X3)
!    if(t+dt/2d0>=tnext) then    !need to load a new file
!      if ( .not. allocated(mlonp)) then    !need to read in the grid data from input file
!        ymdprev=ymd
!        UTsecprev=UTsec
!        ymdnext=ymdprev
!        UTsecnext=UTsecprev
!
!!        if (myid==0) then    !root process
!          !READ IN THE GRID
!          write(filename,*) trim(adjustl(E0dir)),'/simsize.dat'
!          write(*,*) 'Inputting electric field data size from file:  ',trim(adjustl(filename))
!          open(inunit,file=trim(adjustl(filename)),status='old',form='unformatted',access='stream')
!          read(inunit) llon,llat
!          close(inunit)
!          write(*,*) 'Electric field data has llon,llat size:  ',llon,llat
!
!
!!          !MESSAGE WORKERS WITH GRID INFO
!!          do iid=1,lid-1
!!            call mpi_send(llon,1,MPI_INTEGER,iid,tagllon,MPI_COMM_WORLD,ierr)
!!            call mpi_send(llat,1,MPI_INTEGER,iid,tagllat,MPI_COMM_WORLD,ierr)
!!          end do
!          allocate(mlonp(llon),mlatp(llat))    !bit of code duplication with worker code block below...
!
!
!          !IF WE HAVE SINGLETON DIMENSION THEN ALLOCATE SOME SPACE FOR A TEMP
!          !ARRAY FOR INPUTTING INTO INTERP1
!          if (llon==1) then
!            allocate(Edatp(llat))
!          elseif (llat==1) then
!            allocate(Edatp(llon))
!          end if
!
!
!          !NOW READ THE GRID
!          write(filename,*) trim(adjustl(E0dir)),'/simgrid.dat'
!          write(*,*) 'Inputting electric field grid from file:  ',trim(adjustl(filename))
!          open(inunit,file=trim(adjustl(filename)),status='old',form='unformatted',access='stream')
!          read(inunit) mlonp,mlatp
!          close(inunit)
!          write(*,*) 'Electric field data has mlon,mlat extent:  ',minval(mlonp(:)),maxval(mlonp(:)),minval(mlatp(:)), &
!                                                                  maxval(mlatp(:))
!
!!          !NOW SEND THE GRID DATA
!!          do iid=1,lid-1
!!            call mpi_send(mlonp,llon,MPI_DOUBLE_PRECISION,iid,tagmlon,MPI_COMM_WORLD,ierr)    !I think it's okya for these to have the same tag as for the precipitation, but perhaps should be fixed for clarity in the future...
!!            call mpi_send(mlatp,llat,MPI_DOUBLE_PRECISION,iid,tagmlat,MPI_COMM_WORLD,ierr)
!!          end do
!!        else    !workers
!!          call mpi_recv(llon,1,MPI_INTEGER,0,tagllon,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!          call mpi_recv(llat,1,MPI_INTEGER,0,tagllat,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!          allocate(mlonp(llon),mlatp(llat)) 
!!
!!          call mpi_recv(mlonp,llon,MPI_DOUBLE_PRECISION,0,tagmlon,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!          call mpi_recv(mlatp,llat,MPI_DOUBLE_PRECISION,0,tagmlat,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!        end if
!
!        !SPACE TO STORE INPUT DATA
!        allocate(E0xp(llon,llat),E0yp(llon,llat))
!        allocate(E0xiprev(lx2,lx3all),E0xinext(lx2,lx3all),E0yiprev(lx2,lx3all),E0yinext(lx2,lx3all))
!        E0xiprev=0d0; E0yiprev=0d0; E0xinext=0d0; E0yinext=0d0;     !these need to be initialized so that something sensible happens at the beginning
!
!        !ALL PROCESSES NEED TO DEFINED THE OPINTS THAT THEY WILL BE INTERPOLATING ONTO
!        allocate(mloni(lx2*lx3all),mlati(lx2*lx3all))
!        do ix3=1,lx3all
!          do ix2=1,lx2
!            iflat=(ix3-1)*lx2+ix2
!            mlati(iflat)=90d0-x%thetaall(lx1,ix2,ix3)*180d0/pi    !does rool even have this info full grid???
!            mloni(iflat)=x%phiall(lx1,ix2,ix3)*180d0/pi
!          end do
!        end do
!        write(*,*) 'Grid has mlon,mlat range:  ',minval(mloni),maxval(mloni),minval(mlati),maxval(mlati)
!        write(*,*) 'Grid has size:  ',iflat
!      end if
!
!
!      !GRID INFORMATION EXISTS AT THIS POINT SO START READING IN PRECIP DATA
!!      if (myid==0) then    !only root reads file data
!        !read in the data from file
!        write(*,*) 'tprev,tnow,tnext:  ',tprev,t+dt/2d0,tnext
!        ymdtmp=ymdnext
!        UTsectmp=UTsecnext
!        call dateinc(dtE0,ymdtmp,UTsectmp)    !get the date for "next" params
!        filename=date_filename(E0dir,ymdtmp,UTsectmp)     !form the standard data filename
!        write(*,*) 'Pulling electric field data from file:  ',trim(adjustl(filename))
!        open(inunit,file=trim(adjustl(filename)),status='old',form='unformatted',access='stream',iostat=ios)
!        if (ios==0) then    !successful read
!          write(*,*) 'Successfully located input file...'
!          read(inunit) E0xp,E0yp
!        else      !just set everything to zero
!          write(*,*) 'Bad input file, setting everything to some default value...'
!          E0xp=0d0; E0yp=0d0;
!        end if
!        close(inunit)
!
!        write(*,*) 'Min/max values for E0xp:  ',minval(pack(E0xp,.true.)),maxval(pack(E0xp,.true.))
!        write(*,*) 'Min/max values for E0yp:  ',minval(pack(E0yp,.true.)),maxval(pack(E0yp,.true.))
!
!!        !send a full copy of the data to all of the workers
!!        do iid=1,lid-1
!!          write(*,*) 'Sending to worker:  ',iid
!!          call mpi_send(E0xp,llon*llat,MPI_DOUBLE_PRECISION,iid,tagE0xp,MPI_COMM_WORLD,ierr)
!!          call mpi_send(E0yp,llon*llat,MPI_DOUBLE_PRECISION,iid,tagE0yp,MPI_COMM_WORLD,ierr)
!!        end do
!!      else     !workers receive data from root
!!        write(*,*) 'Worker;  ',myid,' receiving...'
!!        call mpi_recv(E0xp,llon*llat,MPI_DOUBLE_PRECISION,0,tagE0xp,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!        call mpi_recv(E0yp,llon*llat,MPI_DOUBLE_PRECISION,0,tagE0yp,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!!      end if
!
!
!      !ALL WORKERS DO SPATIAL INTERPOLATION TO THEIR SPECIFIC GRID SITES
!!      if (myid==0) then
!        write(*,*) 'Initiating precipitation spatial interpolations for date:  ',ymdtmp,' ',UTsectmp
!!      end if
!      if (llon==1) then    !source data has singleton dimension in longitude
!        Edatp=E0xp(1,:)
!        parami=interp1(mlatp,Edatp,mlati)   !will work even for 2D grids, just repeats the data in the lon direction
!        E0xiprev=E0xinext
!        E0xinext=reshape(parami,[lx2,lx3])
!
!        Edatp=E0yp(1,:)
!        parami=interp1(mlatp,Edatp,mlati)
!        E0yiprev=E0yinext
!        E0yinext=reshape(parami,[lx2,lx3])
!      elseif (llat==1) then
!        Edatp=E0xp(:,1)
!        parami=interp1(mlatp,Edatp,mlati)
!        E0xiprev=E0xinext
!        E0xinext=reshape(parami,[lx2,lx3])
!
!        Edatp=E0yp(:,1)
!        parami=interp1(mlatp,Edatp,mlati)
!        E0yiprev=E0yinext
!        E0yinext=reshape(parami,[lx2,lx3])
!      else    !source data is 2D
!        parami=interp2(mlonp,mlatp,E0xp,mloni,mlati)     !interp to temp var.
!        E0xiprev=E0xinext                       !save new pervious
!        E0xinext=reshape(parami,[lx2,lx3all])    !overwrite next with new interpolated input
!
!        parami=interp2(mlonp,mlatp,E0yp,mloni,mlati)
!        E0yiprev=E0yinext
!        E0yinext=reshape(parami,[lx2,lx3all])
!      end if
! !     if (myid==lid/2) then
!        write(*,*) 'Min/max values for E0xi:  ',minval(pack(E0xinext,.true.)),maxval(pack(E0xinext,.true.))
!        write(*,*) 'Min/max values for E0yi:  ',minval(pack(E0yinext,.true.)),maxval(pack(E0yinext,.true.))
!!      end if
!
!
!      !UPDATE OUR CONCEPT OF PREVIOUS AND NEXT TIMES
!      tprev=tnext
!      UTsecprev=UTsecnext
!      ymdprev=ymdnext
!
!      tnext=tprev+dtE0
!      UTsecnext=UTsectmp
!      ymdnext=ymdtmp
!    end if
!
!
!    !INTERPOLATE IN TIME (LINEAR)
!    do ix3=1,lx3all
!      do ix2=1,lx2
!        slope(ix2,ix3)=(E0xinext(ix2,ix3)-E0xiprev(ix2,ix3))/(tnext-tprev)
!        E0xinow(ix2,ix3)=E0xiprev(ix2,ix3)+slope(ix2,ix3)*(t+dt/2d0-tprev)
!
!        slope(ix2,ix3)=(E0yinext(ix2,ix3)-E0yiprev(ix2,ix3))/(tnext-tprev)
!        E0yinow(ix2,ix3)=E0yiprev(ix2,ix3)+slope(ix2,ix3)*(t+dt/2d0-tprev)
!      end do
!    end do
!
!
!    !SOME BASIC DIAGNOSTICS
!!    if (myid==lid/2) then
!      write(*,*) 'tprev,t,tnext:  ',tprev,t+dt/2d0,tnext
!      write(*,*) 'Min/max values for E0xinow:  ',minval(pack(E0xinow,.true.)),maxval(pack(E0xinow,.true.))
!      write(*,*) 'Min/max values for E0yinow:  ',minval(pack(E0yinow,.true.)),maxval(pack(E0yinow,.true.))
!!      write(*,*) 'Min/max values for Qiprev:  ',minval(pack(Qiprev,.true.)),maxval(pack(Qiprev,.true.))
!!      write(*,*) 'Min/max values for E0prev:  ',minval(pack(E0iprev,.true.)),maxval(pack(E0iprev,.true.))
!!      write(*,*) 'Min/max values for Qinext:  ',minval(pack(Qinext,.true.)),maxval(pack(Qinext,.true.))
!!      write(*,*) 'Min/max values for E0next:  ',minval(pack(E0inext,.true.)),maxval(pack(E0inext,.true.))
!!    end if
!
!
!    !ASSIGN RESULTS OF INTERPOLATION TO OUTPUT VARIABLES - THESE NEED TO BE TRANSMITTED BACK TO ROOT...
!!    if (myid==0) then
!!      allocate(E0xinowall(lx2,lx3all),E0yinowall(lx2,lx3all))
!
!!      !Receive a copy of everyone's data
!!      call gather_recv(E0xinow,tagE0xi,E0xinowall)
!!      call gather_recv(E0yinow,tagE0yi,E0yinowall)
!
!      !Place data in the array expected by potential 'solver'
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          do ix1=1,lx1
!            E02all(ix1,ix2,ix3)=E0xinow(ix2,ix3)
!            E03all(ix1,ix2,ix3)=E0yinow(ix2,ix3)      
!          end do
!        end do
!      end do
!
!!      deallocate(E0xinowall,E0yinowall)
!!    else   !workers send their interpolated data to root
!!      call gather_send(E0xinow,tagE0xi)
!!      call gather_send(E0yinow,tagE0yi)
!!    end if
!
!  end subroutine potentialBCs2D_fileinput
!

  subroutine potentialBCs2D(t,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
                                          Vmaxx3,E01all,E02all,E03all,flagdirich)

    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)

    real(wp), intent(in) :: t
    type(curvmesh), intent(in) :: x

    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
    integer, intent(out) :: flagdirich

    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field

    real(wp) :: Phipk
    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
    integer :: im
!    integer, parameter :: lmodes=8
    real(wp) :: phase
    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field

    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt


    !CALCULATE/SET TOP BOUNDARY CONDITIONS
    sigx2=1d0/20d0*(x%x2(lx2)-x%x2(1))
    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
    sigx3=1d0/20d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))

    if (gridflag/=2) then
      Vtopalt=>Vminx1
      Vbotalt=>Vmaxx1
    else
      Vtopalt=>Vmaxx1
      Vbotalt=>Vminx1
    end if

    Phipk=0d0      !pk current density
    flagdirich=0    !Neumann conditions
    do ix3=1,lx3all
      do ix2=1,lx2
        Vtopalt(ix2,ix3)=0d0
      end do
    end do


    !SOME USER INFO
    write(*,*) 'At time:  ',t,'  Max FAC set to be:  ',maxval(pack(abs(Vtopalt),.true.))


    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
    Vbotalt=0d0   !since we need to have no current through bottom boundary
    Vminx2=0d0
    Vmaxx2=0d0
    Vminx3=0d0
    Vmaxx3=0d0


    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
    E01all=0d0
    E02all=0d0
    E03all=0d0

  end subroutine potentialBCs2D


!  subroutine potentialBCs2D_KHI(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!    real(wp), dimension(1:size(E01all,1),1:size(E01all,2),1:size(E01all,3)) ::  vel3    !x3 component of initial drift velocity
!    real(wp), dimension(1:size(E01all,2),1:size(E01all,3)) :: E2slab,Phislab
!
!    real(wp), parameter :: v0=500d0
!    real(wp), parameter :: vn=500d0, voffset=100d0
!    real(wp), parameter :: B1val=-50000d-9    !must match grid structure avg. value
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    meanx3=0d0
!    sigx3=60d3
!    meanx2=0d0
!    sigx2=0.5d3
!    meant=900d0
!    sigt=450d0
!    sigcurv=450d3
!    x30amp=0.1d3    !original successful runs used 2d3, but it was very pronounced...
!    varc=0d0
!    Phipk=120d-3    !pk electric field, 100d-3 works well
!
!
!    !NO CURRENT THROUGH THE TOP BOUNDARY
!    flagdirich=0
!    Vmaxx1=0d0
!
!
!    !CONVERT FAC INTO POTENTIAL DERIVATIVE
!    write(*,*) 'At time:  ',t,'  Max current set to be:  ',maxval(pack(Vmaxx1,.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vminx1=0d0
!    Vminx3=0d0     !these are actually not used if periodic is selected (and it should be for this type of simulation)
!    Vmaxx3=0d0     !also not used if we plan to do a simulation that includes perdiodic boundary conditions which are specified through the input config.dat file
!
!
!    !THE FOLLOWING LINES DUPLICATE CODE FROM THE ICS SUBROUTINE BELOW
!    !FILL IN VELOCITY FIELD (X3 COMP)
!    do ix3=1,lx3all
!      do ix1=1,lx1
!        vel3(ix1,:,ix3)=-v0*tanh(x%x2(1:lx2)/sigx2)+vn+voffset
!      end do
!    end do
!
!
!    !CONVERT TO ELECTRIC FIELD
!    E2slab=vel3(1,:,:)*B1val
!
!
!    !INTEGRATE TO PRODUCE A POTENTIAL OVER GRID
!    Phislab=integral2D1_curv_alt(E2slab,x,1,lx2)
!    do ix1=1,lx1
!      Vmaxx2(ix1,:)=Phislab(lx2,:)
!      Vminx2(ix1,:)=Phislab(1,:)
!    end do
!
!
!   
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!
!  end subroutine potentialBCs2D_KHI
!
!
!  subroutine potentialBCs2D_GDI(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!
!    !SIZES
!!    lx1=size(sig0all,1)
!!    lx2=size(sig0all,2)
!!    lx3all=size(sig0all,3)
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    meanx3=0d0
!    sigx3=60d3
!    meanx2=0d0
!    sigx2=50d3
!    meant=900d0
!    sigt=450d0
!    sigcurv=450d3
!    x30amp=0.1d3    !original successful runs used 2d3, but it was very pronounced...
!    varc=0d0
!
!    Phipk=120d-3    !pk electric field, 100d-3 works well
!
!
!    !NO CURRENT THROUGH THE TOP BOUNDARY
!    flagdirich=0
!    Vmaxx1=0d0
!
!
!    !CONVERT FAC INTO POTENTIAL DERIVATIVE
!    write(*,*) 'At time:  ',t,'  Max current set to be:  ',maxval(pack(Vmaxx1,.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vminx1=0d0
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!    
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=-25e-3
!
!  end subroutine potentialBCs2D_GDI
!
!
!  subroutine potentialBCs2D_3DPCarc2(t,sig0all,x, &
!                   Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !------------------------------------------------------------
!    !-------POPULATES TOP BOUNDARY CONDITIONS FOR POTENTIAL AND
!    !-------SOURCE TERMS FROM BACKGROUND FIELDS, ETC.  THIS 
!    !-------PARTICULAR IMPLEMENTATION USES DIRICHLET CONDITIONS
!    !------------------------------------------------------------
!
!    use phys_consts
!
!    implicit none
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp) :: Jpk
!    integer :: ix1,ix2,ix3,lx1,lx2,lx3all
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc,x2enve    !for setting background field
!
!    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
!
!
!    !SIZES
!    lx1=size(sig0all,1)
!    lx2=size(sig0all,2)
!    lx3all=size(sig0all,3)
!
!
!    !SET POINTERS FOR THIS TYPE OF GRID
!    if (gridflag/=2) then
!      write(*,*) '!!!Using inverted boundary conditions...'
!      Vtopalt=>Vminx1
!      Vbotalt=>Vmaxx1
!    else
!      write(*,*) '!!!Using non-inverted boundary conditions...'
!      Vtopalt=>Vmaxx1
!      Vbotalt=>Vminx1
!    end if
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    flagdirich=0
!
!    sigx2=150d3
!    meanx2=0d0
!    sigx3=25d3
!    meant=120d0
!    sigt=60d0
!    sigcurv=450d3
!    x30amp=0d3
!    varc=200d0
!    Jpk=0.875d-6
!
!
!    if (t<43200d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          meanx3=-varc*meant+varc*t
!          Vtopalt(ix2,ix3)=(Jpk*exp(-(x%x3all(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2)- &
!                          Jpk*exp(-(x%x3all(ix3)-(sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2) )* &
!                              exp(-(x%x2(ix2)-meanx2)**4/2d0/sigx2**4)*exp(-(t-meant)**2/2d0/sigt**2)
!        end do
!      end do
!    else
!      Vtopalt=0d0
!    end if
!
!
!    !CONVERT FAC INTO POTENTIAL DERIVATIVE
!    write(*,*) 'Max FAC set to be:  ',maxval(pack(Vtopalt,.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE EQUIPOTENTIAL WITH TOP EDGE
!    Vbotalt=0d0
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=10d-3
!    E03all=0d0
!
!  end subroutine potentialBCs2D_3DPCarc2
!
!
!  subroutine potentialBCs2D_3DPCarc(t,sig0all,x, &
!                   Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !------------------------------------------------------------
!    !-------POPULATES TOP BOUNDARY CONDITIONS FOR POTENTIAL AND
!    !-------SOURCE TERMS FROM BACKGROUND FIELDS, ETC.  THIS 
!    !-------PARTICULAR IMPLEMENTATION USES DIRICHLET CONDITIONS
!    !------------------------------------------------------------
!
!    use phys_consts
!
!    implicit none
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp) :: Jpk
!    integer :: ix1,ix2,ix3,lx1,lx2,lx3all
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc,x2enve    !for setting background field
!
!    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
!
!
!    !SIZES
!    lx1=size(sig0all,1)
!    lx2=size(sig0all,2)
!    lx3all=size(sig0all,3)
!
!
!    !SET POINTERS FOR THIS TYPE OF GRID
!    if (gridflag/=2) then
!      Vtopalt=>Vminx1
!      Vbotalt=>Vmaxx1
!    else
!      Vtopalt=>Vmaxx1
!      Vbotalt=>Vminx1
!    end if
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    flagdirich=0
!
!    sigx2=150d3
!    meanx2=0d0
!    sigx3=25d3
!    meant=900d0
!    sigt=450d0
!    sigcurv=450d3
!    x30amp=0d3
!    varc=200d0
!    Jpk=0.875d-6
!
!
!    if (t<43200d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          meanx3=-varc*meant+varc*t
!          Vtopalt(ix2,ix3)=(Jpk*exp(-(x%x3all(ix3)-(-sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2)- &
!                          Jpk*exp(-(x%x3all(ix3)-(sigx3+meanx3+x30amp*cos(2d0*pi/(sigx2/2d0)*x%x2(ix2)) ))**2/2d0/sigx3**2) )* &
!                              exp(-(x%x2(ix2)-meanx2)**4/2d0/sigx2**4)*exp(-(t-meant)**2/2d0/sigt**2)
!        end do
!      end do
!    else
!      Vtopalt=0d0
!    end if
!
!
!    !CONVERT FAC INTO POTENTIAL DERIVATIVE
!    write(*,*) 'Max FAC set to be:  ',maxval(pack(Vtopalt,.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE EQUIPOTENTIAL WITH TOP EDGE
!    Vbotalt=0d0
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=10d-3
!    E03all=0d0
!
!  end subroutine potentialBCs2D_3DPCarc
!
!
!  subroutine potentialBCs2D_curved_closed(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    sigx2=1d0/20d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!    sigx3=1d0/20d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!    if (gridflag/=2) then
!      Vtopalt=>Vminx1
!      Vbotalt=>Vmaxx1
!    else
!      Vtopalt=>Vmaxx1
!      Vbotalt=>Vminx1
!    end if
!
!    Phipk=0d0      !pk current density
!    flagdirich=0    !Neumann conditions
!    do ix3=1,lx3all
!      do ix2=1,lx2
!        Vtopalt(ix2,ix3)=0d0
!      end do
!    end do
!
!
!    !SOME USER INFO
!    write(*,*) 'At time:  ',t,'  Max FAC set to be:  ',maxval(pack(abs(Vtopalt),.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vbotalt=0d0   !since we need to have no current through bottom boundary
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!  end subroutine potentialBCs2D_curved_closed
!
!
!  subroutine potentialBCs2D_zeroends3D(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    sigx2=1d0/20d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!    sigx3=1d0/20d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!    if (gridflag/=2) then
!      Vtopalt=>Vminx1
!      Vbotalt=>Vmaxx1
!    else
!      Vtopalt=>Vmaxx1
!      Vbotalt=>Vminx1
!    end if
!
!    Phipk=0d0      !pk current density
!    flagdirich=0    !Neumann conditions
!
!    if (t>1d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vtopalt(ix2,ix3)=Phipk*(exp(-(x%x3all(ix3)-meanx3-sigx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2) - &
!                                 exp(-(x%x3all(ix3)-meanx3+sigx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2))       !note that program always assumes the potential BC for Dirichlet conditions to be stored in Vmaxx1
!        end do
!      end do
!    else
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vtopalt(ix2,ix3)=0d0
!        end do
!      end do
!    end if
!
!
!    !SOME USER INFO
!    write(*,*) 'At time:  ',t,'  Max FAC set to be:  ',maxval(pack(abs(Vtopalt),.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vbotalt=0d0   !since we need to have no current through bottom boundary
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!  end subroutine potentialBCs2D_zeroends3D
!
!
!  subroutine potentialBCs2D_zeropotends(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    sigx2=1d0/10d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!    sigx3=1d0/10d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!    !Phipk=0.20d4    !pk potential, ?????  works well
!    Phipk=0d0
!    flagdirich=1
!
!    if (t>1d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vmaxx1(ix2,ix3)=Phipk*exp(-(x%x3all(ix3)-meanx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2)       !note that program always assumes the potential BC for Dirichlet conditions to be stored in Vmaxx1
!        end do
!      end do
!    else
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vmaxx1(ix2,ix3)=0d0
!        end do
!      end do
!    end if
!
!
!    !SOME USER INFO
!    write(*,*) 'At time:  ',t,'  Max potential set to be:  ',maxval(pack(abs(Vmaxx1),.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vminx1=Vmaxx1   !since we are using Dirichlet conditions, we assume EFL
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!    
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!  end subroutine potentialBCs2D_zeropotends
!
!
!  subroutine potentialBCs2D_Neumann(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out), target :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!    real(wp), dimension(:,:), pointer :: Vtopalt,Vbotalt
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    sigx2=1d0/20d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!    sigx3=1d0/20d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!    if (gridflag/=2) then
!      Vtopalt=>Vminx1
!      Vbotalt=>Vmaxx1
!    else
!      Vtopalt=>Vmaxx1
!      Vbotalt=>Vminx1
!    end if
!
!    Phipk=5d-6      !pk current density
!    flagdirich=0    !Neumann conditions
!
!    if (t>1d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vtopalt(ix2,ix3)=Phipk*(exp(-(x%x3all(ix3)-meanx3-sigx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2) - &
!                                 exp(-(x%x3all(ix3)-meanx3+sigx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2))       !note that program always assumes the potential BC for Dirichlet conditions to be stored in Vmaxx1
!        end do
!      end do
!    else
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vtopalt(ix2,ix3)=0d0
!        end do
!      end do
!    end if
!
!
!    !SOME USER INFO
!    write(*,*) 'At time:  ',t,'  Max FAC set to be:  ',maxval(pack(abs(Vtopalt),.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vbotalt=0d0   !since we need to have no current through bottom boundary
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!  end subroutine potentialBCs2D_Neumann
!
!
!  subroutine potentialBCs2D_Dirichlet(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS IS A SIMPLE GAUSSIAN POTENTIAL PERTURBATION (IN X1,X2,X3 SPAE)
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!    sigx2=1d0/10d0*(x%x2(lx2)-x%x2(1))
!    meanx2=0.5d0*(x%x2(1)+x%x2(lx2))
!    sigx3=1d0/10d0*(x%x3all(lx3all)-x%x3all(1))    !this requires that all workers have a copy of x3all!!!!
!    meanx3=0.5d0*(x%x3all(1)+x%x3all(lx3all))
!
!    Phipk=0.20d4    !pk potential, ?????  works well
!    flagdirich=1
!
!    if (t>1d0) then
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vmaxx1(ix2,ix3)=Phipk*exp(-(x%x3all(ix3)-meanx3)**2/2d0/sigx3**2)* &
!                                 exp(-(x%x2(ix2)-meanx2)**2/2d0/sigx2**2)       !note that program always assumes the potential BC for Dirichlet conditions to be stored in Vmaxx1
!        end do
!      end do
!    else
!      do ix3=1,lx3all
!        do ix2=1,lx2
!          Vmaxx1(ix2,ix3)=0d0
!        end do
!      end do
!    end if
!
!
!    !SOME USER INFO
!    write(*,*) 'At time:  ',t,'  Max potential set to be:  ',maxval(pack(abs(Vmaxx1),.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vminx1=Vmaxx1   !since we are using Dirichlet conditions, we assume EFL
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!    
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    E02all=0d0
!    E03all=0d0
!
!  end subroutine potentialBCs2D_Dirichlet
!
!
!  subroutine potentialBCs2D_GDIcav(t,sig0all,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3, &
!                                          Vmaxx3,E01all,E02all,E03all,flagdirich)
!
!    !THIS WAS USED FOR MY CAVITY-GDI GRL PAPER
!
!    real(wp), intent(in) :: t
!    real(wp), dimension(:,:,:), intent(in) ::  sig0all
!    type(curvmesh), intent(in) :: x
!
!    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
!    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
!    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
!    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
!    integer, intent(out) :: flagdirich
!
!    real(wp), dimension(1:size(Vmaxx1,1),1:size(Vmaxx1,2)) :: Emaxx1    !pseudo-electric field
!
!    real(wp) :: Phipk
!    integer :: ix1,ix2,ix3    !grid sizes are borrow from grid module
!    integer :: im
!    integer, parameter :: lmodes=8
!    real(wp) :: phase
!    real(wp), dimension(1:lx2) :: x3dev    !a little bit surprise we can use grid mod lx2 var as size...
!    real(wp) :: meanx2,sigx2,meanx3,sigx3,meant,sigt,sigcurv,x30amp,varc    !for setting background field
!
!
!    !SIZES
!!    lx1=size(sig0all,1)
!!    lx2=size(sig0all,2)
!!    lx3all=size(sig0all,3)
!
!
!    !CALCULATE/SET TOP BOUNDARY CONDITIONS
!!    meanx3=-320d3    !for 1000km x3 simulation
!!    meanx3=-653d3    !for 1666km x3 simulation
!!    meanx3=-280d3
!!    meanx3=-300d3
!    meanx3=0d0
!!    sigx3=40d3
!!    sigx3=80d3
!!    sigx3=60d3
!    sigx3=60d3
!    meanx2=0d0
!    sigx2=50d3
!    meant=900d0
!    sigt=450d0
!    sigcurv=450d3
!    x30amp=0.1d3    !original successful runs used 2d3, but it was very pronounced...
!    varc=0d0
!
!    Phipk=120d-3    !pk electric field, 100d-3 works well
!
!
!    if (t<90d0) then
!      flagdirich=1
!
!      call random_seed()   !default initialization everytime we are called so we get the same (recomputed) phase each time
!      x3dev=0d0
!      do im=1,lmodes    !comment out if no seed for instability
!        call random_number(phase)
!        phase=phase*2*pi
!        x3dev=x3dev+x30amp*cos(2d0*pi/(4d0*sigx2/real(im,8))*x%x2(1:lx2)+phase)   !center line for patch edge
!      end do
!
!      do ix3=1,lx3all
!        do ix2=1,lx2
!!          Emaxx1(ix2,ix3)=Phipk*exp(-(x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
!!                                exp(-(x2(ix2)-meanx2)**6/2d0/(4d0*sigx2)**6)
!          Emaxx1(ix2,ix3)=Phipk*exp(-(x%x3all(ix3)-meanx3-x3dev(ix2))**10/2d0/sigx3**10)* &
!                                exp(-(x%x2(ix2)-meanx2)**4/2d0/(4d0*sigx2)**4)
!        end do
!      end do
!
!      Vmaxx1=integral2D2_curv_alt(Emaxx1,x,1,lx3all)   !integrate then modulate to avoid edge effects
!      do ix3=1,lx3all
!        do ix2=1,lx2      
!          Vmaxx1(ix2,ix3)=Vmaxx1(ix2,ix3)*exp(-(x%x3all(ix3)-meanx3)**2/2d0/(3d0*sigx3)**2)
!        end do
!      end do
!    else
!      flagdirich=0
!
!      Vmaxx1=0d0
!    end if
!
!
!    !CONVERT FAC INTO POTENTIAL DERIVATIVE
!    write(*,*) 'At time:  ',t,'  Max potential set to be:  ',maxval(pack(Vmaxx1,.true.))
!
!
!    !BOTTOM BOUNDARY IS ALWAYS ZERO CURRENT - SIDES ARE JUST GROUNDED
!    Vminx1=0d0
!    Vminx2=0d0
!    Vmaxx2=0d0
!    Vminx3=0d0
!    Vmaxx3=0d0
!
!    
!    !COMPUTE SOURCE/FORCING TERMS FROM BACKGROUND FIELDS, ETC.
!    E01all=0d0
!    if (t<90d0) then
!      E02all=0d0
!    else
!!      E02all=25d-3
!      E02all=-25d-3
!    end if
!    E03all=0d0
!
!  end subroutine potentialBCs2D_GDIcav
!
end module potentialBCs_mumps
