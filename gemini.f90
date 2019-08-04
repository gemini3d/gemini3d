program gemini

!----------------------------------------------------------
!------THIS IS THE MAIN PROGRAM FOR GEMINI3D
!----------------------------------------------------------

use phys_consts, only : lnchem, lwave, lsp, wp, debug
use grid, only: curvmesh, grid_size,read_grid,clear_grid,lx1,lx2,lx3,lx2all,lx3all
use temporal, only : dt_comm
use timeutils, only: dateinc
use neutral, only : neutral_atmos,make_dneu,neutral_perturb,clear_dneu
use io, only : read_configfile,input_plasma,create_outdir,output_plasma,create_outdir_aur,output_aur
use potential_comm,only : electrodynamics
use multifluid, only : fluid_adv
use mpimod, only : mpisetup, mpibreakdown, mpi_manualgrid, mpigrid, lid, myid
use precipBCs_mod, only: make_precip_fileinput, clear_precip_fileinput
use potentialBCs_mumps, only: clear_potential_fileinput

implicit none

!----------------------------------------------------------
!------VARIABLE DECLARATIONS
!----------------------------------------------------------

!VARIABLES READ IN FROM CONFIG.INI FILE
integer, dimension(3) :: ymd    !year,month,day of simulation
real(wp) :: UTsec      !UT (s)
real(wp) :: UTsec0     !UT start time of simulation (s)
real(wp) :: tdur       !duration of simulation
real(wp), dimension(3) :: activ    !f10.7a,f10.7,ap
real(wp) :: tcfl                       !target CFL number
real(wp) :: Teinf                      !exospheric temperature
integer :: potsolve                   !what type of potential solve
integer :: flagperiodic               !toggles whether or not the grid is treated as periodic in the x3 dimension (affects some of the message passing)
integer :: flagoutput                 !what type of output to do (1 - everything; 2 - avg'd parms.; 3 - ne only)
integer :: flagcap                    !internal capacitance?

!INPUT AND OUTPUT FILES
character(:), allocatable :: infile    !command line argument input file
character(:), allocatable :: outdir    !" " output directory
character(:), allocatable :: indatsize,indatgrid    !grid size and data filenames

!GRID STRUCTURE
type(curvmesh) :: x    !structure containg grid locations, finite differences, etc.:  see grid module for details

!STATE VARIABLES
real(wp), dimension(:,:,:,:), allocatable :: ns,vs1,vs2,vs3,Ts    !fluid state variables
real(wp), dimension(:,:,:), allocatable :: E1,E2,E3,J1,J2,J3      !electrodynamic state variables
real(wp), dimension(:,:,:), allocatable :: rhov2,rhov3,B1,B2,B3   !inductive state vars. (for future use - except for B1 which is used for the background field)
real(wp), dimension(:,:,:), allocatable :: rhom,v1,v2,v3          !inductive auxiliary
real(wp), dimension(:,:,:,:), allocatable :: nn                   !neutral density array
real(wp), dimension(:,:,:), allocatable :: Tn,vn1,vn2,vn3         !neutral temperature and velocities
real(wp), dimension(:,:,:), allocatable :: Phiall                 !full-grid potential solution.  To store previous time step value
real(wp), dimension(:,:,:), allocatable :: iver                   !integrated volume emission rate of aurora calculated by GLOW

!TEMPORAL VARIABLES
real(wp) :: t=0._wp,dt=1e-6_wp,dtprev      !time from beginning of simulation (s) and time step (s)
real(wp) :: tout,dtout    !time for next output and time between outputs
real(wp) :: tstart,tfin   !temp. vars. for measuring performance of code blocks
integer :: it,isp        !time and species loop indices

!WORK ARRAYS
real(wp), allocatable :: dl1,dl2,dl3     !these are grid distances in [m] used to compute Courant numbers

!NEUTRAL PERTURBATION VARIABLES
integer :: flagdneu                  !toggles neutral perturbations (0 - none; 1 - file-based neutral inputs)
integer :: interptype                !toggles whether the neutral input data are interpreted (0 - Cartesian; 1 - axisymmetric)
real(wp) :: drhon,dzn                 !finite differences for the neutral input data in the horizontal and vertical directions
real(wp) :: sourcemlat,sourcemlon     !mag. lat./long for the neutral source location
character(:), allocatable :: sourcedir          !directory where neutral input data are located
real(wp) :: dtneu                     !time interval [s] in between neutral inputs

!PRECIPITATION FILE INPUT VARIABLES
integer :: flagprecfile              ! flag toggling precipitation file input (0 - no; 1 - yes)
real(wp) :: dtprec                    ! time interval between precip. inputs
character(:), allocatable :: precdir ! directory containing precip. input files

!ELECTRIC FIELD FILE INPUT VARIABLES
integer :: flagE0file                ! flag toggling electric field (potential BCs) file input (0 - no; 1 - yes)
real(wp) :: dtE0                      ! time interval between electric field file inputs
character(:), allocatable :: E0dir   ! directory containing electric field file input data

!GLOW MODULE INPUT VARIABLES
integer :: flagglow                     !flag toggling GLOW module run (include aurora) (0 - no; 1 - yes)
real(wp) :: dtglow                      !time interval between GLOW runs (s)
real(wp) :: dtglowout                   !time interval between GLOW auroral outputs (s)
real(wp) :: tglowout                    !time for next GLOW output

!FOR HANDLING OUTPUT
integer :: argc
character(256) :: argv
integer :: lid2in,lid3in


!TO CONTROL THROTTLING OF TIME STEP
real(wp), parameter :: dtscale=2d0

!----------------------------------------------------------
!------MAIN PROGRAM
!----------------------------------------------------------

argc = command_argument_count()
if (argc < 2) error stop 'must specify .ini file to configure simulation and output directory'

!INITIALIZE MESSING PASSING VARIABLES, IDS ETC.
call mpisetup()
print *, 'Process:  ',myid,' of:  ',lid-1,' online...'


!READ FILE INPUT
call get_command_argument(1,argv)
infile = trim(argv)

call read_configfile(infile, ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf,potsolve,flagperiodic,flagoutput,flagcap, &
                     indatsize,indatgrid,flagdneu,interptype,sourcemlat,sourcemlon,dtneu,drhon,dzn,sourcedir,flagprecfile, &
                     dtprec,precdir,flagE0file,dtE0,E0dir,flagglow,dtglow,dtglowout)

!!CHECK THE GRID SIZE AND ESTABLISH A PROCESS GRID
call grid_size(indatsize)

select case (argc)
  case (4,5) !user specified process grid
    call get_command_argument(3,argv)
    read(argv,*) lid2in
    call get_command_argument(4,argv)
    read(argv,*) lid3in
    call mpi_manualgrid(lx2all,lx3all,lid2in,lid3in)
    if (argc == 5) then
      call get_command_argument(5,argv)
      if (argv == '-d' .or. argv == '-debug')  debug = .true.
    endif
  case default
    if (argc == 3) then
      call get_command_argument(3,argv)
      if (argv == '-d' .or. argv == '-debug')  debug = .true.
    endif
    call mpigrid(lx2all,lx3all)    !following grid_size these are in scope
end select


!LOAD UP THE GRID STRUCTURE/MODULE VARS. FOR THIS SIMULATION
call read_grid(indatsize,indatgrid,flagperiodic,x)     !read in a previously generated grid from filenames listed in input file


!CREATE/PREP OUTPUT DIRECTORY AND OUTPUT SIMULATION SIZE AND GRID DATA; ONLY THE ROOT PROCESS WRITES OUTPUT DATA
call get_command_argument(2,argv)
outdir = trim(argv)

if (myid==0) then
  call create_outdir(outdir,infile,indatsize,indatgrid,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
  if (flagglow/=0) call create_outdir_aur(outdir)
end if


!ALLOCATE ARRAYS (AT THIS POINT ALL SIZES ARE SET FOR EACH PROCESS SUBGRID)
allocate(ns(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs1(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp),vs2(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), &
  vs3(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp), Ts(-1:lx1+2,-1:lx2+2,-1:lx3+2,lsp))
allocate(rhov2(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhov3(-1:lx1+2,-1:lx2+2,-1:lx3+2),B1(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
         B2(-1:lx1+2,-1:lx2+2,-1:lx3+2),B3(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(v1(-1:lx1+2,-1:lx2+2,-1:lx3+2),v2(-1:lx1+2,-1:lx2+2,-1:lx3+2), &
         v3(-1:lx1+2,-1:lx2+2,-1:lx3+2),rhom(-1:lx1+2,-1:lx2+2,-1:lx3+2))
allocate(E1(lx1,lx2,lx3),E2(lx1,lx2,lx3),E3(lx1,lx2,lx3),J1(lx1,lx2,lx3),J2(lx1,lx2,lx3),J3(lx1,lx2,lx3))
allocate(nn(lx1,lx2,lx3,lnchem),Tn(lx1,lx2,lx3),vn1(lx1,lx2,lx3), vn2(lx1,lx2,lx3),vn3(lx1,lx2,lx3))
call make_dneu()    !allocate space for neutral perturbations in case they are used with this run
call make_precip_fileinput()


!ALLOCATE MEMORY FOR ROOT TO STORE CERTAIN VARS. OVER ENTIRE GRID
if (myid==0) then
  allocate(Phiall(lx1,lx2all,lx3all))
end if

!ALLOCATE MEMORY FOR AURORAL EMISSIONS, IF CALCULATED
if (flagglow/=0) then
  allocate(iver(lx2,lx3,lwave))
end if

!LOAD ICS AND DISTRIBUTE TO WORKERS (REQUIRES GRAVITY FOR INITIAL GUESSING)
call input_plasma(x%x1,x%x2all,x%x3all,indatsize,ns,vs1,Ts)


!ROOT/WORKERS WILL ASSUME THAT THE MAGNETIC FIELDS AND PERP FLOWS START AT ZERO
!THIS KEEPS US FROM HAVING TO HAVE FULL-GRID ARRAYS FOR THESE STATE VARS (EXCEPT
!FOR IN OUTPUT FNS.).  IF A SIMULATIONS IS DONE WITH INTERTIAL CAPACITANCE THERE
!WILL BE A FINITE AMOUNT OF TIME FOR THE FLOWS TO 'START UP', BUT THIS SHOULDN'T
!BE TOO MUCH OF AN ISSUE.  WE ALSO NEED TO SET THE BACKGROUND MAGNETIC FIELD STATE
!VARIABLE HERE TO WHATEVER IS SPECIFIED IN THE GRID STRUCTURE (THESE MUST BE CONSISTENT)
rhov2=0d0; rhov3=0d0; v2=0d0; v3=0d0;
B2=0d0; B3=0d0;
B1(1:lx1,1:lx2,1:lx3)=x%Bmag      !this assumes that the grid is defined s.t. the x1 direction corresponds to the magnetic field direction (hence zero B2 and B3).


!INITIALIZE ELECTRODYNAMIC QUANTITIES FOR POLARIZATION CURRENT
if (myid==0) then
  Phiall=0d0     !only root store entire potential array
end if
E1=0d0; E2=0d0; E3=0d0;
vs2=0d0; vs3=0d0;

!INITIALIZE AURORAL EMISSION MAP
if(flagglow/=0) then
  iver=0.0_wp
end if

!MAIN LOOP
UTsec=UTsec0; it=1; t=0d0; tout=t; tglowout=t;
do while (t<tdur)
  !TIME STEP CALCULATION
  dtprev=dt
  call dt_comm(t,tout,tglowout,flagglow,tcfl,ns,Ts,vs1,vs2,vs3,B1,B2,B3,x,potsolve,dt)
  if (it>1) then
    if(dt/dtprev > dtscale) then     !throttle how quickly we allow dt to increase
      dt=dtscale*dtprev
      if (myid==0) then
        print *, 'Throttling dt to:  ',dt
      end if
    end if
  end if


  !COMPUTE BACKGROUND NEUTRAL ATMOSPHERE USING MSIS00.  PRESENTLY THIS ONLY GETS CALLED
  !ON THE FIRST TIME STEP DUE TO A NEED TO KEEP A CONSTANT BACKGROUND (I.E. ONE NOT VARYING
  !IN TIME) FOR SOME SIMULATIONS USING NEUTRAL INPUT.  REALISTICALLY THIS NEEDS TO BE
  !RECALLED EVERY SO OFTEN (MAYBE EVERY 10-15 MINS)
  if (it==1) then
    call cpu_time(tstart)
    call neutral_atmos(ymd,UTsec,x%glat,x%glon,x%alt,activ,nn,Tn)
    vn1=0d0; vn2=0d0; vn3=0d0     !hard-code these to zero for the first time step
    call cpu_time(tfin)

    if (myid==0) then
      print *, 'Neutral background calculated in time:  ',tfin-tstart
    end if
  end if


  !GET NEUTRAL PERTURBATIONS FROM ANOTHER MODEL
  if (flagdneu==1) then
    call cpu_time(tstart)
    if (it==1) then    !this triggers the code to load the neutral frame correspdonding ot the beginning time of the simulation
      if (myid==0) print *, '!!!Attempting initial load of neutral dynamics files!!!' // &
                              ' This is a workaround that fixes the restart code...',t-dt
      call neutral_perturb(interptype,dt,dtneu,t-dtneu,ymd,UTsec-dtneu,sourcedir,drhon,dzn, &
                                  sourcemlat,sourcemlon,x,nn,Tn,vn1,vn2,vn3)
    end if
    call neutral_perturb(interptype,dt,dtneu,t,ymd,UTsec,sourcedir,drhon,dzn,sourcemlat,sourcemlon,x,nn,Tn,vn1,vn2,vn3)
    call cpu_time(tfin)
    if (myid==0 .and. debug) print *, 'Neutral perturbations calculated in time:  ',tfin-tstart
  end if

  !POTENTIAL SOLUTION
  call cpu_time(tstart)
  call electrodynamics(it,t,dt,nn,vn2,vn3,Tn,sourcemlat,ns,Ts,vs1,B1,vs2,vs3,x, &
                        potsolve,flagcap,E1,E2,E3,J1,J2,J3, &
                        Phiall,flagE0file,dtE0,E0dir,ymd,UTsec)
  call cpu_time(tfin)
  if (myid==0 .and. debug) print *, 'Electrodynamics total solve time:  ',tfin-tstart


  !UPDATE THE FLUID VARIABLES
  call cpu_time(tstart)
  call fluid_adv(ns,vs1,Ts,vs2,vs3,J1,E1,Teinf,t,dt,x,nn,vn1,vn2,vn3,Tn,iver,activ(2),activ(1),ymd,UTsec, &
                 flagprecfile,dtprec,precdir,flagglow,dtglow)
  call cpu_time(tfin)
  if (myid==0 .and. debug) print *, 'Multifluid total solve time:  ',tfin-tstart


  !NOW OUR SOLUTION IS FULLY UPDATED SO UPDATE TIME VARIABLES TO MATCH...
  it=it+1; t=t+dt;
  if (myid==0 .and. debug) print *, 'Moving on to time step (in sec):  ',t,'; end time of simulation:  ',tdur
  call dateinc(dt,ymd,UTsec)
  if (myid==0) print *, 'Current date',ymd,'Current UT time:  ',UTsec


  !OUTPUT
  if (abs(t-tout) < 1d-5) then   !close enough to warrant an output now...
    call cpu_time(tstart)
    call output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
    call cpu_time(tfin)
    if (myid==0 .and. debug) print *, 'Plasma output done for time step:  ',t,' in cpu_time of:  ',tfin-tstart

    tout=tout+dtout
  end if

  !GLOW OUTPUT
  if ((flagglow/=0).and.(abs(t-tglowout) < 1d-5)) then !same as plasma output
    call cpu_time(tstart)
    call output_aur(outdir,flagglow,ymd,UTsec,iver)
    call cpu_time(tfin)
    if (myid==0) then
      print *, 'Auroral output done for time step:  ',t,' in cpu_time of: ',tfin-tstart
    end if
    tglowout=tglowout+dtglowout
  end if
end do


!DEALLOCATE MAIN PROGRAM DATA
deallocate(ns,vs1,vs2,vs3,Ts)
deallocate(E1,E2,E3,J1,J2,J3)
deallocate(nn,Tn,vn1,vn2,vn3)

if (myid==0) then
  deallocate(Phiall)
end if

if (flagglow/=0) then
  deallocate(iver)
end if

!DEALLOCATE MODULE VARIABLES (MAY HAPPEN AUTOMATICALLY IN F2003???)
call clear_grid(x)
call clear_dneu()
call clear_precip_fileinput()
call clear_potential_fileinput()


!SHUT DOWN MPI
call mpibreakdown()

end program gemini
