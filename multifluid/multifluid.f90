module multifluid

use phys_consts, only : gammas,kB,ms,mindensdiv,mindens,mindensnull
use mpimod
use grid
use ionization
use sources
use calculus
use advec_mpi
use diffusion
use precipBCs_mod
use temporal, only : sza
implicit none

integer, parameter :: lprec=2    !number of precipitating electron populations


contains

  subroutine fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,Teinf,t,dt,x,nn,vn1,vn2,vn3,Tn,f107,f107a,ymd,UTsec, &
                       flagprecfile,dtprec,precdir)    !J1 needed for heat conduction; E1 for momentum equation

    !------------------------------------------------------------
    !-------THIS SUBROUTINE ADVANCES ALL OF THE FLUID VARIABLES 
    !------ BY TIME STEP DT.
    !------------------------------------------------------------ 

    real(8), dimension(-1:,-1:,-1:,:), intent(inout) ::  ns,vs1,Ts
    real(8), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3
    real(8), dimension(:,:,:), intent(in) :: J1       !needed for thermal conduction in electron population
    real(8), dimension(:,:,:), intent(inout) :: E1    !will have ambipolar field added into it in this procedure...

    real(8), intent(in) :: Teinf,t,dt

    type(curvmesh), intent(in) :: x                   !grid structure variable

    real(8), dimension(:,:,:,:), intent(in) :: nn
    real(8), dimension(:,:,:), intent(in) :: vn1,vn2,vn3,Tn
    real(8), intent(in) :: f107,f107a
    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec

    integer, intent(in) :: flagprecfile
    real(8), intent(in) :: dtprec
    character(*), intent(in) :: precdir

    integer :: isp
    real(8) :: tstart,tfin

    real(8), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2,size(ns,4)) ::  rhovs1,rhoes
    real(8), dimension(-1:size(ns,1)-2,-1:size(ns,2)-2,-1:size(ns,3)-2) :: param
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: A,B,C,D,E,paramtrim,rhoeshalf,lambda,beta,chrgflux
    real(8), dimension(0:size(ns,1)-3,0:size(ns,2)-3,0:size(ns,3)-3) :: divvs
    real(8), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4) :: v1i
    real(8), dimension(1:size(vs1,1)-4,1:size(vs1,2)-3,1:size(vs1,3)-4) :: v2i
    real(8), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-3) :: v3i

    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)) :: Pr,Lo
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)-1) :: Prprecip,Prpreciptmp
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: Qeprecip
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: chi
    real(8), dimension(1:size(ns,2)-4,1:size(ns,3)-4,lprec) :: W0,PhiWmWm2

    integer :: iprec
    real(8), dimension(1:size(vs1,1)-3,1:size(vs1,2)-4,1:size(vs1,3)-4) :: v1iupdate    !temp interface velocities for art. viscosity
    real(8), dimension(1:size(vs1,1)-4,1:size(vs1,2)-4,1:size(vs1,3)-4) :: dv1iupdate    !interface diffs. for art. visc.
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,size(ns,4)) :: Q
    real(8), parameter :: xicon=3d0    !decent value for closed field-line grids extending to high altitudes.  


    !CALCULATE THE INTERNAL ENERGY AND MOMENTUM FLUX DENSITIES (ADVECTION AND SOURCE SOLUTIONS ARE DONE IN THESE VARIABLES)
    do isp=1,lsp
      rhovs1(:,:,:,isp)=ns(:,:,:,isp)*ms(isp)*vs1(:,:,:,isp)
      rhoes(:,:,:,isp)=ns(:,:,:,isp)*kB*Ts(:,:,:,isp)/(gammas(isp)-1d0)
    end do


    !ADVECTION SUBSTEP (CONSERVED VARIABLES SHOULD BE UPDATED BEFORE ENTERING)
    call cpu_time(tstart)
    chrgflux=0d0
    do isp=1,lsp
      call advec_prep_mpi(isp,x%flagper,ns,rhovs1,vs1,vs2,vs3,rhoes,v1i,v2i,v3i)    !role-agnostic communcation pattern (all-to-neighbors)

      if(isp<lsp) then   !electron info found from charge neutrality and current density
        param=ns(:,:,:,isp)
        param=advec3D_MC_mpi(param,v1i,v2i,v3i,dt,x,0)   !last argument is tensor rank of thing being advected
        ns(:,:,:,isp)=param
  
        param=rhovs1(:,:,:,isp)
        param=advec3D_MC_mpi(param,v1i,v2i,v3i,dt,x,1)
        rhovs1(:,:,:,isp)=param
  
        vs1(:,:,:,isp)=rhovs1(:,:,:,isp)/(ms(isp)*max(ns(:,:,:,isp),mindensdiv))
        chrgflux=chrgflux+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*vs1(1:lx1,1:lx2,1:lx3,isp)
      else
        ns(:,:,:,lsp)=sum(ns(:,:,:,1:lsp-1),4)
  !      vs1(1:lx1,1:lx2,1:lx3,lsp)=1d0/ns(1:lx1,1:lx2,1:lx3,lsp)/qs(lsp)*(J1-chrgflux)   !density floor needed???
        vs1(1:lx1,1:lx2,1:lx3,lsp)=-1d0/max(ns(1:lx1,1:lx2,1:lx3,lsp),mindensdiv)/qs(lsp)*chrgflux   !really not strictly correct, should include current density
      end if
  
      param=rhoes(:,:,:,isp)
      param=advec3D_MC_mpi(param,v1i,v2i,v3i,dt,x,0)
      rhoes(:,:,:,isp)=param
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Completed advection substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if


    !CLEAN DENSITY AND VELOCITY - SETS THE NULL CELLS TO SOME SENSIBLE VALUE SO
    !THEY DON'T MESS UP FINITE DIFFERENCES LATER
    call clean_param(x,1,ns)
    call clean_param(x,2,vs1)
 

    !ARTIFICIAL VISCOSITY (NOT REALLY NEED BELOW 1000 KM ALT.).  NOTE THAT WE DON'T CHECK WHERE SUBCYCLING IS NEEDED SINCE, IN MY EXPERIENCE THEN CODE IS BOMBING ANYTIME IT IS...
    do isp=1,lsp-1
      v1iupdate(1:lx1+1,:,:)=0.5d0*(vs1(0:lx1,1:lx2,1:lx3,isp)+vs1(1:lx1+1,1:lx2,1:lx3,isp))    !compute an updated interface velocity (only in x1-direction)
      dv1iupdate=v1iupdate(2:lx1+1,:,:)-v1iupdate(1:lx1,:,:)
      Q(:,:,:,isp)=ns(1:lx1,1:lx2,1:lx3,isp)*ms(isp)*0.25d0*xicon**2*(min(dv1iupdate,0d0))**2   !note that viscosity does not have/need ghost cells
    end do
    Q(:,:,:,lsp)=0d0
 
  
    !NONSTIFF/NONBALANCE INTERNAL ENERGY SOURCES (RK2 INTEGRATION) 
    call cpu_time(tstart)
    do isp=1,lsp
      call RK2_prep_mpi(isp,x%flagper,vs1,vs2,vs3)    !role-agnostic mpi, all-to-neighbor
  
      divvs=div3D(vs1(0:lx1+1,0:lx2+1,0:lx3+1,isp),vs2(0:lx1+1,0:lx2+1,0:lx3+1,isp), &
                 vs3(0:lx1+1,0:lx2+1,0:lx3+1,isp),x,0,lx1+1,0,lx2+1,0,lx3+1)    !diff with one set of ghost cells to preserve second order accuracy over the grid
      paramtrim=rhoes(1:lx1,1:lx2,1:lx3,isp)
      rhoeshalf=paramtrim-dt/2d0*(paramtrim*(gammas(isp)-1d0)+Q(:,:,:,isp))*divvs(1:lx1,1:lx2,1:lx3) !t+dt/2 value of internal energy, use only interior points of divvs for second order accuracy
  
      paramtrim=paramtrim-dt*(rhoeshalf*(gammas(isp)-1d0)+Q(:,:,:,isp))*divvs(1:lx1,1:lx2,1:lx3)
      rhoes(1:lx1,1:lx2,1:lx3,isp)=paramtrim
  
      Ts(:,:,:,isp)=(gammas(isp)-1d0)/kB*rhoes(:,:,:,isp)/max(ns(:,:,:,isp),mindensdiv)
      Ts(:,:,:,isp)=max(Ts(:,:,:,isp),100d0)
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Completed compression substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
  

    !CLEAN TEMPERATURE
    call clean_param(x,3,Ts)


    !DIFFUSION OF ENERGY
    call cpu_time(tstart)
    do isp=1,lsp
      param=Ts(:,:,:,isp)     !temperature for this species
      call thermal_conduct(isp,param,ns(:,:,:,isp),nn,J1,lambda,beta)
  
      call diffusion_prep(isp,x,lambda,beta,ns(:,:,:,isp),param,A,B,C,D,E,Tn,Teinf)
!      param=backEuler3D(param,A,B,C,D,E,dt,dx1,dx1i)    !1st order method, likely deprecated but needs to be kept here for debug purposes, perhaps?
      param=TRBDF23D(param,A,B,C,D,E,dt,x)
      Ts(:,:,:,isp)=param
      Ts(:,:,:,isp)=max(Ts(:,:,:,isp),100d0) 
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Completed energy diffusion substep for time step:  ',t,' in cpu_time of:  ',tfin-tstart
    end if
 

    !ZZZ - CLEAN TEMPERATURE BEFORE CONVERTING TO INTERNAL ENERGY
    call clean_param(x,3,Ts)
    do isp=1,lsp
      rhoes(:,:,:,isp)=ns(:,:,:,isp)*kB*Ts(:,:,:,isp)/(gammas(isp)-1d0)
    end do


    !LOAD ELECTRON PRECIPITATION PATTERN
    if (flagprecfile==1) then
      call precipBCs_fileinput(dt,dtprec,t,ymd,UTsec,precdir,x,W0,PhiWmWm2)
    else     !no file input specified, so just call 'regular' function
      call precipBCs(t,x,W0,PhiWmWm2)
    end if
  
    !STIFF/BALANCED ENERGY SOURCES
    call cpu_time(tstart)
    Prprecip=0d0
    if (gridflag/=0) then
      do iprec=1,lprec    !loop over the different populations of precipitation (2 here?), accumulating production rates
        Prpreciptmp=ionrate_fang08(W0(:,:,iprec),PhiWmWm2(:,:,iprec),x%alt,nn,Tn)    !calculation based on Fang et al [2008]
        Prprecip=Prprecip+Prpreciptmp
      end do
    else    !do not compute impact ionization on a closed mesh (presumably there is no source of energetic electrons at these lats.)
      if (myid==0) then
        write(*,*) 'Looks like we have a closed grid, so skipping impact ionization for time step:  ',t
      end if
    end if

    !now add in photoionization sources
    chi=sza(ymd,UTsec,x%glat,x%glon)
    if (myid==0) then
      write(*,*) 'Computing photoionization for time:  ',t,' using sza range of (root only):  ', &
                  minval(chi)*180d0/pi,maxval(chi)*180d0/pi
    end if
    Prpreciptmp=photoionization(x,nn,Tn,chi,f107,f107a)
    if (myid==0) then
      write(*,*) 'Min/max root production rates for time:  ',t,' :  ',minval(pack(Prpreciptmp,.true.)), &
                  maxval(pack(Prpreciptmp,.true.))
    end if
    Prprecip=Prprecip+Prpreciptmp
    Prprecip=max(Prprecip,1d-5)    !enforce some minimum production rate to preserve conditioning for species that rely on constant production, some testing should probably be done to see what the best choice is...

    Qeprecip=eheating(nn,Tn,Prprecip,ns)   !thermal electron heating rate from Swartz and Nisbet, (1978)
  
    call srcsEnergy(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo)
    do isp=1,lsp
      if (isp==lsp) then
        Pr(:,:,:,lsp)=Pr(:,:,:,lsp)+Qeprecip
      end if
      paramtrim=rhoes(1:lx1,1:lx2,1:lx3,isp)
      paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
      rhoes(1:lx1,1:lx2,1:lx3,isp)=paramtrim
  
      Ts(:,:,:,isp)=(gammas(isp)-1d0)/kB*rhoes(:,:,:,isp)/max(ns(:,:,:,isp),mindensdiv)
      Ts(:,:,:,isp)=max(Ts(:,:,:,isp),100d0)
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Energy sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
    end if


    !CLEAN TEMPERATURE
    call clean_param(x,3,Ts)

  
    !ALL VELOCITY SOURCES
    call cpu_time(tstart)
    call srcsMomentum(nn,vn1,Tn,ns,vs1,vs2,vs3,Ts,E1,Q,x,Pr,Lo)    !added artificial viscosity...
    do isp=1,lsp-1
      paramtrim=rhovs1(1:lx1,1:lx2,1:lx3,isp)
      paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
      rhovs1(1:lx1,1:lx2,1:lx3,isp)=paramtrim
  
      vs1(:,:,:,isp)=rhovs1(:,:,:,isp)/(ms(isp)*max(ns(:,:,:,isp),mindensdiv))
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Velocity sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
    end if


    !ELECTRON VELOCITY SOLUTION
    chrgflux=0d0
    do isp=1,lsp-1
      chrgflux=chrgflux+ns(1:lx1,1:lx2,1:lx3,isp)*qs(isp)*vs1(1:lx1,1:lx2,1:lx3,isp)
    end do
  !  vs1(1:lx1,1:lx2,1:lx3,lsp)=1d0/max(ns(1:lx1,1:lx2,1:lx3,lsp),mindensdiv)/qs(lsp)*(J1-chrgflux)   !density floor needed???
    vs1(1:lx1,1:lx2,1:lx3,lsp)=-1d0/max(ns(1:lx1,1:lx2,1:lx3,lsp),mindensdiv)/qs(lsp)*chrgflux    !don't bother with FAC contribution...


    !CLEAN VELOCITY
    call clean_param(x,2,vs1)

 
    !ALL MASS SOURCES
    call cpu_time(tstart)
    call srcsContinuity(nn,vn1,vn2,vn3,Tn,ns,vs1,vs2,vs3,Ts,Pr,Lo)
    Pr(:,:,:,1:6)=Pr(:,:,:,1:6)+Prprecip
    do isp=1,lsp-1
      paramtrim=ns(1:lx1,1:lx2,1:lx3,isp)
      paramtrim=ETD_uncoupled(paramtrim,Pr(:,:,:,isp),Lo(:,:,:,isp),dt)
      ns(1:lx1,1:lx2,1:lx3,isp)=paramtrim    !should there be a density floor here???  I think so...
    end do
    call cpu_time(tfin)
    if (myid==0) then
      write(*,*) 'Mass sources substep for time step:  ',t,'done in cpu_time of:  ',tfin-tstart
    end if
 
 
    !ELECTRON DENSITY SOLUTION
    ns(:,:,:,lsp)=sum(ns(:,:,:,1:lsp-1),4)
  

    !CLEAN DENSITY (CONSERVED VARIABLES WILL BE RECOMPUTED AT THE BEGINNING OF NEXT TIME STEP
    call clean_param(x,1,ns)

    !should the electron velocity be recomputed here now that densities have changed...
  
  end subroutine fluid_adv


  subroutine clean_param(x,paramflag,param)

    !------------------------------------------------------------
    !-------THIS SUBROUTINE ZEROS OUT ALL NULL CELLS AND HANDLES
    !-------POSSIBLE NULL ARTIFACTS AT BOUNDARIES
    !------------------------------------------------------------ 

    type(curvmesh), intent(in) :: x
    integer, intent(in) :: paramflag
    real(8), dimension(-1:,-1:,-1:,:), intent(inout) :: param     !note that this is 4D and is meant to include ghost cells

    real(8), dimension(-1:size(param,1)-2,-1:size(param,2)-2,-1:size(param,3)-2,lsp) :: paramnew
    integer :: isp,ix1,ix2,ix3,iinull,ix1beg,ix1end


    select case (paramflag)
      case (1)    !density
        param(:,:,:,1:lsp-1)=max(param(:,:,:,1:lsp-1),mindens)
        param(:,:,:,lsp)=sum(param(:,:,:,1:lsp-1),4)       !enforce charge neutrality based on ion densities

        do isp=1,lsp             !set null cells to some value
          if (isp==1) then
            do iinull=1,x%lnull
              ix1=x%inull(iinull,1)
              ix2=x%inull(iinull,2)
              ix3=x%inull(iinull,3)

              param(ix1,ix2,ix3,isp)=mindensnull*1d-2
            end do
          else
            do iinull=1,x%lnull
              ix1=x%inull(iinull,1)
              ix2=x%inull(iinull,2)
              ix3=x%inull(iinull,3)
          
              param(ix1,ix2,ix3,isp)=mindensnull
            end do
          end if
        end do
      case (2)    !velocity
        do isp=1,lsp       !set null cells to zero mometnum
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)

            param(ix1,ix2,ix3,isp)=0d0
          end do
        end do

        !FORCE THE BORDER CELLS TO BE SAME AS THE FIRST INTERIOR CELL (deals with some issues on dipole grids)
        do isp=1,lsp
          do ix3=1,lx3
            do ix2=1,lx2
              ix1beg=1
              do while(x%nullpts(ix1beg,ix2,ix3)<0.5d0 .and. ix1beg<lx1)     !find the first non-null index for this field line, need to be careful if no null points exist...
                ix1beg=ix1beg+1
              end do
          
              ix1end=ix1beg
              do while(x%nullpts(ix1end,ix2,ix3)>0.5d0 .and. ix1end<lx1)     !find the first non-null index for this field line
                ix1end=ix1end+1
              end do
 
              if (ix1beg /= lx1) then    !only do this if we actually have null grid points 
                param(ix1beg,ix2,ix3,isp)=param(ix1beg+1,ix2,ix3,isp)
                param(ix1end,ix2,ix3,isp)=param(ix1end-1,ix2,ix3,isp)
              end if
            end do
          end do
        end do
      case (3)    !temperature
        param=max(param,100d0)     !temperature floor

        do isp=1,lsp       !set null cells to some value
          do iinull=1,x%lnull
            ix1=x%inull(iinull,1)
            ix2=x%inull(iinull,2)
            ix3=x%inull(iinull,3)
      
            param(ix1,ix2,ix3,isp)=100d0
          end do
        end do 
      case default    !do nothing...
        write(*,*)  '!non-standard parameter selected in clean_params...'
    end select

  end subroutine clean_param

end module multifluid
