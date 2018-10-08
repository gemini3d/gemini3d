program magcalc

!----------------------------------------------------------
!------THIS IS THE MAIN PROGRAM FOR COMPUTING MAGNETIC FIELDS
!------FROM OUTPUT FROM A SIMULATIONS DONE BY GEMINI3D.  THIS
!------PROGRAM VERY MUCH MIRRORS THE SETUP OF THE MAIN GEMINI.F90
!------CODE.  
!----------------------------------------------------------

use phys_consts, only : lnchem,pi,mu0
use grid
use temporal, only : dateinc
use io, only : read_configfile,input_plasma_currents,create_outdir_mag,output_magfields
use mpimod

implicit none

!----------------------------------------------------------
!------VARIABLE DECLARATIONS
!----------------------------------------------------------

!VARIABLES READ IN FROM CONFIG.DAT FILE
integer, dimension(3) :: ymd    !year,month,day of simulation
real(8) :: UTsec      !UT (s)
real(8) :: UTsec0     !UT start time of simulation (s)
real(8) :: tdur       !duration of simulation
real(8), dimension(3) :: activ    !f10.7a,f10.7,ap
real(8) :: tcfl                       !target CFL number
real(8) :: Teinf                      !exospheric temperature
integer :: potsolve                   !what type of potential solve 
integer :: flagperiodic               !toggles whether or not the grid is treated as periodic in the x3 dimension (affects some of the message passing)
integer :: flagoutput                 !what type of output to do (1 - everything; 2 - avg'd parms.; 3 - ne only)
integer :: flagcap                    !internal capacitance?

!INPUT AND OUTPUT FILES
character(:), allocatable :: infile    !command line argument input file
character(:), allocatable :: outdir    !" " output directory
character(:), allocatable :: indatsize,indatgrid    !grid size and data filenames
character(:), allocatable :: fieldpointfile
integer :: u

!GRID STRUCTURE
type(curvmesh) :: x    !structure containg grid locations, finite differences, etc.:  see grid module for details

!STATE VARIABLES
real(8), dimension(:,:,:), allocatable :: J1,J2,J3      !electrodynamic state variables

!TEMPORAL VARIABLES
real(8) :: t=0d0,dt      !time from beginning of simulation (s) and time step (s)
real(8) :: tout,dtout    !time for next output and time between outputs
real(8) :: tstart,tfin   !temp. vars. for measuring performance of code blocks
integer :: it,isp        !time and species loop indices

!WORK ARRAYS
integer :: flag2D
real(8), dimension(:,:,:), allocatable :: xp,yp,zp     !source coordinates computed from simulation sub-grid
real(8), dimension(:), allocatable  :: xf,yf,zf        !field point coordinates (flat list)
real(8), dimension(:), allocatable :: r,theta,phi
real(8), dimension(:,:,:), allocatable :: dV
integer :: ipoints,lpoints    !size of field point arrays
real(8), dimension(:,:,:), allocatable :: proj_e1er,proj_e2er,proj_e3er
real(8), dimension(:,:,:), allocatable :: proj_e1etheta,proj_e2etheta,proj_e3etheta
real(8), dimension(:,:,:), allocatable :: proj_e1ephi,proj_e2ephi,proj_e3ephi
real(8), dimension(:,:,:), allocatable :: Jx,Jy,Jz
real(8), dimension(:,:,:), allocatable :: Rx,Ry,Rz,Rcubed
real(8), dimension(:,:,:), allocatable :: integrand,integrandavg
real(8), dimension(:,:,:), allocatable :: alt
real(8), dimension(:), allocatable :: Br,Btheta,Bphi
real(8), dimension(:), allocatable :: Brall,Bthetaall,Bphiall
integer :: ix1,ix2,ix3
real(8) :: rmean,thetamean

!NEUTRAL PERTURBATION VARIABLES (UNUSED)
integer :: flagdneu                  !toggles neutral perturbations (0 - none; 1 - file-based neutral inputs)
integer :: interptype                !toggles whether the neutral input data are interpreted (0 - Cartesian; 1 - axisymmetric)
real(8) :: drhon,dzn                 !finite differences for the neutral input data in the horizontal and vertical directions
real(8) :: sourcemlat,sourcemlon     !mag. lat./long for the neutral source location
character(:), allocatable :: sourcedir          !directory where neutral input data are located
real(8) :: dtneu                     !time interval [s] in between neutral inputs

!PRECIPITATION FILE INPUT VARIABLES (UNUSED)
integer :: flagprecfile              ! flag toggling precipitation file input (0 - no; 1 - yes)
real(8) :: dtprec                    ! time interval between precip. inputs
character(:), allocatable :: precdir ! directory containing precip. input files

!ELECTRIC FIELD FILE INPUT VARIABLES (UNUSED)
integer :: flagE0file                ! flag toggling electric field (potential BCs) file input (0 - no; 1 - yes)
real(8) :: dtE0                      ! time interval between electric field file inputs
character(:), allocatable :: E0dir   ! directory containing electric field file input data

!FOR HANDLING INPUT
integer :: argc
character(256) :: argv

!----------------------------------------------------------
!------MAIN PROGRAM
!----------------------------------------------------------

argc = command_argument_count()
if (argc < 2) error stop 'magcalc.f90 --> must specify .ini file to configure simulation and location of &
                          output data from simulation (i.e. plasma parameters and fields) and &
                          and input field point file'


!INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()
write(*,*) 'magcalc.f90 --> Process:  ',myid,' of:  ',lid-1,' online...'


!READ CONFIG FILE FROM OUTPUT DIRECTORY
call get_command_argument(1,argv)
outdir=trim(argv)
infile=outdir//'/inputs/config.ini'
!infile = trim(argv)
if (myid==0) then
  write(*,*) 'Simulation data directory:  ',outdir
  write(*,*) 'Input config file:  ',infile
end if
call read_configfile(infile,ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf,potsolve,flagperiodic,flagoutput,flagcap, &
                     indatsize,indatgrid,flagdneu,interptype,sourcemlat,sourcemlon,dtneu,drhon,dzn,sourcedir,flagprecfile, &
                     dtprec,precdir,flagE0file,dtE0,E0dir)


!LOAD UP THE GRID STRUCTURE/MODULE VARS. FOR THIS SIMULATION - THIS ALSO PERMUTES DIMENSIONS OF 2D GRID, IF NEEDED
call read_grid(indatsize,indatgrid,flagperiodic,x)     !read in a previously generated grid from filename listed in input file, distribute subgrids to individual workers
if (lx2==1) then
  flag2D=1    !a 2D simulation was done, which changes how the integrations go...
else
  flag2D=0
end if


!GRAB THE INFO FOR WHERE THE OUTPUT CALCULATIONS ARE STORED
!call get_command_argument(2,argv)
!outdir = trim(argv)    !this should be the base output directory that the simulation results have been stored in
call get_command_argument(2,argv)
fieldpointfile=trim(argv)    !this file contains the field points at which we are computing magnetic perturbations, it will be copied into the output directory


!SET UP DIRECTORY TO STORE OUTPUT FILES
if (myid==0) then
  call create_outdir_mag(outdir,fieldpointfile)
end if


!ALLOCATE ARRAYS (AT THIS POINT ALL SIZES ARE SET FOR EACH PROCESS SUBGRID)
allocate(J1(lx1,lx2,lx3),J2(lx1,lx2,lx3),J3(lx1,lx2,lx3))


!NOW DEAL WITH THE UNMPRIMED COORDINATES
open(newunit=u,file=fieldpointfile,form='unformatted',access='stream',action='read')
read(u) lpoints    !size of coordinates for field points
if (myid==0) write(*,*) 'magcalc.f90 --> Number of field points:  ',lpoints
allocate(r(lpoints),theta(lpoints),phi(lpoints))
read(u) r,theta,phi
if (myid==0) write(*,*) 'magcalc.f90 --> Range of r,theta,phi',minval(r),maxval(r),minval(theta), &
                           maxval(theta),minval(phi),maxval(phi)
rmean=sum(r)/size(r)
thetamean=sum(theta)/size(theta)
allocate(xf(lpoints),yf(lpoints),zf(lpoints))
xf(:)=r(:)
!yf(:)=r(:)*theta(:)
!zf(:)=r(:)*sin(theta(:))*phi(:)
yf(:)=rmean*theta(:)
zf(:)=rmean*sin(thetamean)*phi(:)


!GET POSITIONS (CARTESIAN) SET UP FOR MAGNETIC COMPUTATIONS.  THESE ARE PRIMED COORDINATES (SOURCE COORDS, I.E. THE SIM GRID)
allocate(xp(lx1,lx2,lx3),yp(lx1,lx2,lx3),zp(lx1,lx2,lx3))
xp(:,:,:)=x%alt(:,:,:)+Re                               !radial distance from Earth's center
!yp(:,:,:)=xp(:,:,:)*x%theta(:,:,:)                      !southward distance (in the direction of the theta spherical coordinate)
!zp(:,:,:)=xp(:,:,:)*sin(x%theta(:,:,:))*x%phi(:,:,:)    !eastward distance
yp(:,:,:)=rmean*x%theta(:,:,:)
zp(:,:,:)=rmean*sin(thetamean)*x%phi(:,:,:)


!COMPUTE A SOURCE DIFFERENTIAL VOLUME FOR INTEGRALS
allocate(dV(lx1,lx2,lx3))
if (flag2D/=1) then   !3D differential volume
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        dV(ix1,ix2,ix3)=x%h1(ix1,ix2,ix3)*x%h2(ix1,ix2,ix3)*x%h3(ix1,ix2,ix3)*x%dx1(ix1)*x%dx2(ix2)*x%dx3(ix3)
      end do
    end do
  end do
else         !plane geometry assumption
  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        dV(ix1,ix2,ix3)=x%h1(ix1,ix2,ix3)*x%h3(ix1,ix2,ix3)*x%dx1(ix1)*x%dx3(ix3)
      end do
    end do
  end do
end if


!COMPUTE NEEDED PROJECTIONS
allocate(proj_e1er(lx1,lx2,lx3),proj_e2er(lx1,lx2,lx3),proj_e3er(lx1,lx2,lx3))
allocate(proj_e1etheta(lx1,lx2,lx3),proj_e2etheta(lx1,lx2,lx3),proj_e3etheta(lx1,lx2,lx3))
allocate(proj_e1ephi(lx1,lx2,lx3),proj_e2ephi(lx1,lx2,lx3),proj_e3ephi(lx1,lx2,lx3))
allocate(alt(lx1,lx2,lx3))
alt(:,:,:)=x%alt
proj_e1er(:,:,:)=sum(x%e1*x%er,4)
proj_e2er(:,:,:)=sum(x%e2*x%er,4)
proj_e3er(:,:,:)=sum(x%e3*x%er,4)
proj_e1etheta(:,:,:)=sum(x%e1*x%etheta,4)
proj_e2etheta(:,:,:)=sum(x%e2*x%etheta,4)
proj_e3etheta(:,:,:)=sum(x%e3*x%etheta,4)
proj_e1ephi(:,:,:)=sum(x%e1*x%ephi,4)
proj_e2ephi(:,:,:)=sum(x%e2*x%ephi,4)
proj_e3ephi(:,:,:)=sum(x%e3*x%ephi,4)


!DEALLOCATE GRID MODULE VARIABLES TO SAVE MEMORY (PROGRAM DOESN'T ACTUALLY NEED THESE ONCE X,Y,Z CREATED
call clear_grid(x)
deallocate(r,theta,phi)


!STORAGE FOR MAGNETIC FIELD CALCULATIONS
allocate(Rx(lx1,lx2,lx3),Ry(lx1,lx2,lx3),Rz(lx1,lx2,lx3))
allocate(Rcubed(lx1,lx2,lx3))
allocate(integrand(lx1,lx2,lx3),integrandavg(lx1-1,max(lx2-1,1),lx3-1))    !latter is cell centered hence -1 in size, max is needed to prevent zero sized array
allocate(Br(lpoints),Btheta(lpoints),Bphi(lpoints))
allocate(Brall(lpoints),Bthetaall(lpoints),Bphiall(lpoints))    !only used by root, but I think workers need to have space allocated for this


!MAIN LOOP
UTsec=UTsec0; it=1; t=0d0; tout=t;
do while (t<tdur)
  !TIME STEP CALCULATION
  dt=dtout    !only compute magnetic field at times when we've done output


  !READ IN THE FULL PLASMA AND FIELD DATA FROM THE OUTPUT FILE (NOTE THAT WE NEED TO KNOW OUTPUT TYPE DONE)
  if (it==1) then
    UTsec=UTsec+0.000001d0
  end if
  if (it==2) then
    UTsec=UTsec-0.000001d0
  end if
  call input_plasma_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)    !now everyone has their piece of data


  !FORCE PARALLEL CURRENTS TO ZERO BELOW 80KM
  if(myid==0) write(*,*) 'Nullifying low altitude currents...'
  where (alt<75d3)
    J1=0d0
  end where


  !DEAL WITH THE WEIRD EDGE ARTIFACTS THAT WE GET IN THE PARALLEL CURRENT
  !SOMETIMES
  if(myid==0) write(*,*) 'Fixing potential edge artifacts...'
  if (myid==lid-1) then
    if (lx3>2) then    !do a ZOH
      J1(:,:,lx3-1)=J1(:,:,lx3-2)
      J1(:,:,lx3)=J1(:,:,lx3-2)
    else
      J1(:,:,lx3-1)=0d0
      J1(:,:,lx3)=0d0
    end if
  end if
  if (myid==0) then
    if (lx3>2) then    !do a ZOH
      J1(:,:,1)=J1(:,:,3)
      J1(:,:,2)=J1(:,:,3)
    else
      J1(:,:,1)=0d0
      J1(:,:,2)=0d0
    end if
  end if


  !ROTATE MAGNETIC FIELDS INTO VERTICAL,SOUTH,EAST COMPONENTS
  if (myid==0) then
    write(*,*) 'magcalc.f90 --> Rotating currents into geomagnetic coordinates...'
  end if
  Jx=J1*proj_e1er+J2*proj_e2er+J3*proj_e3er                 !vertical
  Jy=J1*proj_e1etheta+J2*proj_e2etheta+J3*proj_e3etheta     !south
  Jz=J1*proj_e1ephi+J2*proj_e2ephi+J3*proj_e3ephi           !east
!  write(*,*) ' Currents min/max:  ',maxval(abs(Jx)),maxval(abs(Jy)),maxval(abs(Jz))


  !COMPUTE MAGNETIC FIELDS
  do ipoints=1,lpoints
    if (myid == 0) then
      write(*,*) 'magcalc.f90 --> Computing magnetic field for field point:  ',ipoints,' out of:  ',lpoints
      write(*,*) '            --> ...for time:  ',ymd,UTsec
    end if


    Rx(:,:,:)=xf(ipoints)-xp(:,:,:)
    Ry(:,:,:)=yf(ipoints)-yp(:,:,:)
    Rz(:,:,:)=zf(ipoints)-zp(:,:,:)
!    write(*,*)  'min/max distances:  ',maxval(abs(Rx)),maxval(abs(Ry)),maxval(abs(Rz))
!    write(*,*) mu0,pi

    if (flag2D/=1) then
      Rcubed(:,:,:)=(Rx**2+Ry**2+Rz**2)**(3d0/2d0)   !this really is R**3

      !Bx calculation
      integrand(:,:,:)=mu0/4d0/pi*(Jy*Rz-Jz*Ry)/Rcubed
      integrandavg(:,:,:)=1d0/8d0*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )
      Br(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))    !without dim= input it just sums everything which is what we want

      !By
      integrand(:,:,:)=-mu0/4d0/pi*(Jx*Rz-Jz*Rx)/Rcubed
      integrandavg(:,:,:)=1d0/8d0*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )
      Btheta(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))    !without dim= input it just sums everything which is what we want


      !Bz
      integrand(:,:,:)=mu0/4d0/pi*(Jx*Ry-Jy*Rx)/Rcubed
      integrandavg(:,:,:)=1d0/8d0*( integrand(1:lx1-1,1:lx2-1,1:lx3-1) + integrand(2:lx1,1:lx2-1,1:lx3-1) + &
                           integrand(1:lx1-1,2:lx2,1:lx3-1) + integrand(2:lx1,2:lx2,1:lx3-1) + &
                           integrand(1:lx1-1,1:lx2-1,2:lx3) + integrand(2:lx1,1:lx2-1,2:lx3) + &
                           integrand(1:lx1-1,2:lx2,2:lx3) + integrand(2:lx1,2:lx2,2:lx3) )
      Bphi(ipoints)=sum(integrandavg*dV(2:lx1,2:lx2,2:lx3))    !without dim= input it just sums everything which is what we want
    else
      Rcubed(:,:,:)=Rx**2+Ry**2    !not really R**3 in 2D, just the denominator of the integrand

      !Bx
      integrand(:,:,:)=mu0/4d0/pi*(-2d0*Jz*Ry)/Rcubed
      integrandavg(:,:,:)=1d0/4d0*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx2) + integrand(2:lx1,:,2:lx3) )
      Br(ipoints)=sum(integrandavg*dV(2:lx1,:,2:lx3))    !without dim= input it just sums everything which is what we want

      !By
      integrand(:,:,:)=mu0/4/pi*(2*Jz*Rx)/Rcubed
      integrandavg(:,:,:)=1d0/4d0*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx2) + integrand(2:lx1,:,2:lx3) )
      Btheta(ipoints)=sum(integrandavg*dV(2:lx1,:,2:lx3))    !without dim= input it just sums everything which is what we want

      !Bz
      integrand(:,:,:)=mu0/4d0/pi*2d0*(Jx*Ry-Jy*Rx)/Rcubed
      integrandavg(:,:,:)=1d0/4d0*( integrand(1:lx1-1,:,1:lx3-1) + integrand(2:lx1,:,1:lx3-1) + &
                           integrand(1:lx1-1,:,2:lx2) + integrand(2:lx1,:,2:lx3) )
      Bphi(ipoints)=sum(integrandavg*dV(2:lx1,:,2:lx3))    !without dim= input it just sums everything which is what we want
    end if
  end do

!  write(*,*) '  --> Min/max values of field',minval(Br),maxval(Br),minval(Btheta),maxval(Btheta), &
!                                             minval(Bphi),maxval(Bphi)


  !A REDUCE OPERATION IS NEEDED HERE TO COMBINE MAGNETIC FIELDS (LINEAR SUPERPOSITION) FROM ALL WORKERS
  if (myid ==0) then
    write(*,*) 'Attempting reduction of magnetic field...'
  end if
  call mpi_reduce(Br,Brall,lpoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(Btheta,Bthetaall,lpoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(Bphi,Bphiall,lpoints,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (myid == 0) then
    write(*,*) 'magcalc.f90 --> Reduced magnetic field...'
    write(*,*) '  --> Min/max values of reduced field',minval(Brall),maxval(Brall),minval(Bthetaall),maxval(Bthetaall), &
                                               minval(Bphiall),maxval(Bphiall)
  end if


  !OUTPUT SHOULD BE DONE FOR EVERY INPUT FILE THAT HAS BEEN READ IN
  if (myid==0) then
    call cpu_time(tstart)
    call output_magfields(outdir,ymd,UTsec,Brall,Bthetaall,Bphiall)   !mag field data already reduced so just root needs to output
    call cpu_time(tfin)
    write(*,*) 'magcalc.f90 --> Output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart
  end if
  tout=tout+dtout


  !NOW OUR SOLUTION IS FULLY UPDATED SO UPDATE TIME VARIABLES TO MATCH...
  it=it+1; t=t+dt;
  if (myid==0) then
    write(*,*) 'magcalc.f90 --> Moving on to time step (in sec):  ',t,'; end time of simulation:  ',tdur
  end if
  call dateinc(dt,ymd,UTsec)
  if (myid==0) then
    write(*,*) 'magcalc.f90 --> Current date',ymd,'Current UT time:  ',UTsec
  end if
end do


!DEALLOCATE MAIN PROGRAM DATA
deallocate(J1,J2,J3)
deallocate(xp,yp,zp)
deallocate(xf,yf,zf)
deallocate(dV)
deallocate(Rx,Ry,Rz)
deallocate(proj_e1er,proj_e2er,proj_e3er)
deallocate(proj_e1etheta,proj_e2etheta,proj_e3etheta)
deallocate(proj_e1ephi,proj_e2ephi,proj_e3ephi)
deallocate(Rcubed)
deallocate(integrand,integrandavg)
deallocate(Br,Btheta,Bphi)


!SHUT DOWN MPI
call mpibreakdown()

end program magcalc
