module temporal

!DO NOT FIX THESE WARNINGS - THEY ARE FOR UNUSED VARIABLES THAT MAY BE LEVERAGED
!IN LATER RELEASES
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused
!parameter ‘ns’ [-Wunused-parameter]
!   pure subroutine
!dt_calc(tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,dx1i,dx2i,dx3i,potsolve,cour1,cour2,cour3,dt) ^
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b1’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b2’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘b3’ [-Wunused-parameter]
!/home/zettergm/zettergmdata/GEMINI/temporal/temporal.f90:65:0: warning: unused parameter ‘potsolve’ [-Wunused-parameter]

use phys_consts, only:  kB,mu0,ms,lsp,pi
use mpimod
use grid, only:  curvmesh

implicit none

contains


  subroutine dt_comm(t,tout,tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,potsolve,dt)

    real(wp), intent(in) :: t,tout,tcfl
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1,vs2,vs3
    real(wp),  dimension(-1:,-1:,-1:), intent(in) :: B1,B2,B3
    type(curvmesh), intent(in) :: x
    integer, intent(in) :: potsolve
    real(wp), intent(out) :: dt

    real(wp), dimension(lsp) :: cour1,cour2,cour3
    integer :: iid,isp
    real(wp) :: dttmp


    call dt_calc(tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x%dl1i,x%dl2i,x%dl3i,potsolve,cour1,cour2,cour3,dt)
    if (myid/=0) then
      call mpi_send(dt,1,MPI_DOUBLE_PRECISION,0,tagdt,MPI_COMM_WORLD,ierr)   !send what I think dt should be
      call mpi_recv(dt,1,MPI_DOUBLE_PRECISION,0,tagdt,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)   !receive roots decision   
    else
      !FIGURE OUT GLOBAL DT REQUIRED FOR STABILITY
      do iid=1,lid-1
        call mpi_recv(dttmp,1,MPI_DOUBLE_PRECISION,iid,tagdt,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  
        if (dttmp < dt) then
          dt=dttmp
        end if
      end do
  
      !CHECK WHETHER WE'D OVERSTEP OUR TARGET OUTPUT TIME
      if (t+dt>tout) then
        dt=tout-t
      end if
  
      !DON'T ALLOW ZERO DT
      dt=max(dt,1d-6)
  
      !SEND GLOBAL DT TO ALL WORKERS
      do iid=1,lid-1
        call mpi_send(dt,1,MPI_DOUBLE_PRECISION,iid,tagdt,MPI_COMM_WORLD,ierr)
      end do
  
      write(*,*) 'dt figured to be:  ',dt
      write(*,*) 'x1,x2,x3 courant numbers (root process only!):  '
      do isp=1,lsp
        write(*,'(a4,f4.2,a2,f4.2,a2,f4.2)') '    ',cour1(isp),', ',cour2(isp),', ',cour3(isp)    !these are roots courant numbers
      end do
      write(*,*) 'Min and max density:  ',minval(pack(ns(:,:,:,7),.true.)),maxval(pack(ns(:,:,:,7),.true.))
    end if
  
  end subroutine dt_comm


  pure subroutine dt_calc(tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,dx1i,dx2i,dx3i,potsolve,cour1,cour2,cour3,dt)

    !------------------------------------------------------------
    !-------COMPUTE TIME STEP SUCH THAT STABILITY CONDITION IS
    !-------SATISFIED.  NOTE THAT THE DIFFERENTIALS ARE ASSUMED
    !-------TO HAVE UNITS OF DISTANCE, SO THEY MUST IMPLICITLY
    !-------INCLUDE THE METRIC FACTORS.
    !------------------------------------------------------------


    real(wp), intent(in) :: tcfl
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1,vs2,vs3
    real(wp),  dimension(-1:,-1:,-1:), intent(in) :: B1,B2,B3
    real(wp), dimension(:,:,:), intent(in) :: dx1i
    real(wp), dimension(:,:,:), intent(in) :: dx2i
    real(wp), dimension(:,:,:), intent(in) :: dx3i
    integer, intent(in) :: potsolve
    real(wp), dimension(lsp), intent(out) :: cour1,cour2,cour3
    real(wp), intent(out) :: dt

    real(wp), dimension(lsp) :: gridrate1,gridrate2,gridrate3
    real(wp) :: vsnd
    real(wp) :: rhom,Bmag,vA
    integer :: lx1,lx2,lx3,ix1,ix2,ix3,isp

    lx1=size(Ts,1)-4
    lx2=size(Ts,2)-4
    lx3=size(Ts,3)-4

    gridrate1=0d0
    gridrate2=0d0
    gridrate3=0d0


    !EVALUATE TIME STEP AGAINST LOCAL SOUND SPEED AND ADVECTION
    do isp=1,lsp
      do ix3=1,lx3
        do ix2=1,lx2
          do ix1=1,lx1
            if (isp<lsp) then
              vsnd=sqrt(kB*Ts(ix1,ix2,ix3,isp)/ms(isp)+5d0/3d0*kB*Ts(ix1,ix2,ix3,lsp)/ms(isp))
            else
              vsnd=0d0
            end if

            gridrate1(isp)=max((vsnd+abs(vs1(ix1,ix2,ix3,isp)))/dx1i(ix1,ix2,ix3),gridrate1(isp))
            gridrate2(isp)=max(abs(vs2(ix1,ix2,ix3,isp))/dx2i(ix1,ix2,ix3),gridrate2(isp))
            gridrate3(isp)=max(abs(vs3(ix1,ix2,ix3,isp))/dx3i(ix1,ix2,ix3),gridrate3(isp))
          end do
        end do
      end do
    end do


!    !CHECK GRIDRATE MAX AGAINST LOCAL ALFVEN SPEED (IF THIS SIMULATION IS INDUCTIVE)
!    !NOTE THAT THIS SHOULD REALLY INCLUDE MAGNETOSONIC (FAST) MODES IN GRIDRATE DETERMINATION
!    !AS OF 2/10/2016 THIS CODE IS NOT USED AT ALL, BUT IS KEPT FOR FUTURE DEVELOPMENT
!    if (potsolve == 2) then
!      do ix3=1,lx3
!        do ix2=1,lx2
!          do ix1=1,lx1
!            rhom=0d0
!            do isp=1,lsp
!              rhom=rhom+ms(isp)*ns(ix1,ix2,ix3,isp)
!            end do
!            Bmag=sqrt(B1(ix1,ix2,ix3)**2+B2(ix1,ix2,ix3)**2+B3(ix1,ix2,ix3)**2)
!
!            vA=Bmag/sqrt(rhom*mu0)
!
!            do isp=1,lsp
!              gridrate1(isp)=max(vA/dx1i(ix1),gridrate1(isp))
!              gridrate2(isp)=max(vA/dx2i(ix2),gridrate2(isp))
!              gridrate3(isp)=max(vA/dx3i(ix3),gridrate3(isp))
!            end do
!          end do
!        end do
!      end do
!    end if


    !ENFORCE A MINIMUM VALUE FOR THE GRIDRATE (WHICH IS CFL/DT, GRID POINTS PER SECOND)
    gridrate1=max(gridrate1,1d-10)
    gridrate2=max(gridrate2,1d-10)
    gridrate3=max(gridrate3,1d-10)

    dt=tcfl*min(minval(1d0/gridrate1),minval(1d0/gridrate2),minval(1d0/gridrate3))
    cour1=dt*gridrate1
    cour2=dt*gridrate2
    cour3=dt*gridrate3
  end subroutine dt_calc


  subroutine dateinc(dt,ymd,UTsec)

    real(wp), intent(in) :: dt
    integer, dimension(3), intent(inout) :: ymd
    real(wp), intent(inout) :: UTsec

    integer :: year,month,day
    integer :: monthinc          !let's us know whether we incremented the month
    integer :: daymax     !number of days February should have for given year

    year=ymd(1); month=ymd(2); day=ymd(3);


    UTsec=UTsec+dt
    if (UTsec>=86400d0) then
      UTsec=mod(UTsec,86400d0)
      day=day+1          !roll the day over

      select case(month)     !check whether the month needs to be rolled over
        case(4,6,9,11)    !months having 30 days
          daymax=30
        case(2)    !Feb.
          if (mod(year,4)==0) then     !leap year
            daymax=29
          else          !normal year
            daymax=28
          end if 
        case default    !eveyone else has 31 days
          daymax=30
      end select

      if (day>daymax) then    !roll the month over
        day=1
        month=month+1
        monthinc=1
      else
        monthinc=0
      end if

      if (monthinc==1 .and. month>12) then    !roll the year over
        month=1
        year=year+1
      end if
    end if

    ymd(1)=year; ymd(2)=month; ymd(3)=day;    !replace input array with new date

  end subroutine dateinc


  pure integer function doy_calc(ymd)

    integer, dimension(3), intent(in) :: ymd

    integer, dimension(12) :: monthday


    monthday=[31,28,31,20,31,30,31,31,30,31,30,31]
    if (mod(ymd(1),4)==0) then
      monthday(2)=29     !leap year adjustment
    end if
    if(ymd(2)==1) then
      doy_calc=ymd(3)
    else
      doy_calc=sum(monthday(1:ymd(2)-1))+ymd(3) 
    end if

  end function doy_calc


  function sza(ymd,UTsec,glat,glon)   !computes sza in radians

    !------------------------------------------------------------
    !-------CALCULATE SOLAR ZENITH ANGLE OVER A GIVEN GET OF LAT/LON
    !------------------------------------------------------------

    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:), intent(in) :: glat,glon

    real(wp), dimension(size(glat,1),size(glat,2),size(glat,3)) :: sza

    real(wp) :: doy,soldecrad
    real(wp), dimension(size(glat,1),size(glat,2),size(glat,3)) :: lonrad,LThrs,latrad,hrang


    !SOLAR DECLINATION ANGLE
    doy=doy_calc(ymd)
    soldecrad=-23.44d0*cos(2d0*pi/365d0*(doy+10d0))*pi/180d0;


    !HOUR ANGLE
    lonrad=glon*pi/180d0
    where (lonrad>pi)
      lonrad=lonrad-2d0*pi
    end where
    where (lonrad<-2d0*pi)
      lonrad=lonrad+2d0*pi
    end where
    LThrs=UTsec/3600d0+lonrad/(pi/12d0)
    hrang=(12d0-LThrs)*(pi/12d0)


    !SOLAR ZENITH ANGLE
    latrad=glat*pi/180d0
    sza=acos(sin(soldecrad)*sin(latrad)+cos(soldecrad)*cos(latrad)*cos(hrang))

  end function sza

end module temporal

