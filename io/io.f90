module io

!HANDLES INPUT AND OUTPUT OF PLASMA STATE PARAMETERS (NOT GRID INPUTS)
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only : kB,ms,pi,lsp, wp
use calculus
use mpimod
use grid, only : gridflag,flagswap,lx1,lx2,lx3,lx3all
implicit none


!NONE OF THESE VARIABLES SHOULD BE ACCESSED BY PROCEDURES OUTSIDE THIS MODULE
character(:), allocatable, private :: indatfile                    !initial condition data files from input configuration file


contains


  subroutine read_configfile(infile,ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf,potsolve,flagperiodic, &
                     flagoutput,flagcap,indatsize,indatgrid,flagdneu,interptype, &
                     sourcemlat,sourcemlon,dtneu,drhon,dzn,sourcedir,flagprecfile,dtprec,precdir, &
                     flagE0file,dtE0,E0dir)

    !------------------------------------------------------------
    !-------READS THE INPUT CONFIGURAITON FILE ANDE ASSIGNS
    !-------VARIABLES FOR FILENAMES, SIZES, ETC.
    !------------------------------------------------------------

    character(*), intent(in) :: infile
    integer, dimension(3), intent(out):: ymd
    real(8), intent(out) :: UTsec0
    real(8), intent(out) :: tdur
    real(8), intent(out) :: dtout
    real(8), dimension(3), intent(out) :: activ
    real(8), intent(out) :: tcfl
    real(8), intent(out) :: Teinf
    integer, intent(out) :: potsolve, flagperiodic, flagoutput, flagcap 
    integer, intent(out) :: flagdneu
    integer, intent(out) :: interptype
    real(8), intent(out) :: sourcemlat,sourcemlon
    real(8), intent(out) :: dtneu
    real(8), intent(out) :: drhon,dzn
    integer, intent(out) :: flagprecfile
    real(8), intent(out) :: dtprec
    character(:), allocatable, intent(out) :: indatsize,indatgrid, sourcedir, precdir, E0dir
    integer, intent(out) :: flagE0file
    real(8), intent(out) :: dtE0


    character(256) :: buf
    integer :: u


    !READ CONFIG.DAT FILE FOR THIS SIMULATION
    open(newunit=u, file=infile, status='old', action='read')
    read(u,*) ymd(3),ymd(2),ymd(1)
    read(u,*) UTsec0
    read(u,*) tdur
    read(u,*) dtout
    read(u,*) activ(1),activ(2),activ(3)
    read(u,*) tcfl
    read(u,*) Teinf
    read(u,*) potsolve
    read(u,*) flagperiodic
    read(u,*) flagoutput
    read(u,*) flagcap
    read(u,'(a256)') buf  !apparently format specifier needed, else it reads just one character
    indatsize = trim(buf)
    read(u,'(a256)') buf
    indatgrid = trim(buf)
    read(u,'(a256)') buf
    indatfile = trim(buf)

    !PRINT SOME DIAGNOSIC INFO FROM ROOT
    if (myid==0) then
      print '(A,I6,A1,I0.2,A1,I0.2)', infile//': simulation ymd is:  ',ymd(1),'/',ymd(2),'/',ymd(3)
      print '(A51,F10.3)', 'start time is:  ',UTsec0
      print '(A51,F10.3)', 'duration is:  ',tdur
      print '(A51,F10.3)', 'output every:  ',dtout
      print *, '...using input data files:  '
      print *, '  ',indatsize
      print *, '  ',indatgrid
      print *, '  ',indatfile
    end if


    !NEUTRAL PERTURBATION INPUT INFORMATION
    read(u,*) flagdneu
    if( flagdneu==1) then
      read(u,*) interptype
      read(u,*) sourcemlat,sourcemlon
      read(u,*) dtneu
      read(u,*) drhon,dzn
      read(u,'(A256)') buf
      sourcedir=trim(adjustl(buf))
      if (myid ==0) then
        print *, 'Neutral disturbance mlat,mlon:  ',sourcemlat,sourcemlon
        print *, 'Neutral disturbance cadence (s):  ',dtneu
        print *, 'Neutral grid resolution (m):  ',drhon,dzn
        print *, 'Neutral disturbance data files located in directory:  ',sourcedir
      end if
    else                              !just set it to something
      interptype=0
      sourcemlat=0._wp; sourcemlon=0._wp;
      dtneu=0._wp
      drhon=0._wp; dzn=0._wp;
      sourcedir=''
    end if

    !PRECIPITATION FILE INPUT INFORMATION
    read(u,*) flagprecfile
    if (flagprecfile==1) then    !get the location of the precipitation input files
      read(u,*) dtprec

      read(u,'(A256)') buf
      precdir = trim(adjustl(buf))  ! NOTE: Need adjustl() since trim() is only for TRAILING spaces
      
      if (myid==0) then
        print '(A,F10.3)', 'Precipitation file input cadence (s):  ',dtprec
        print *, 'Precipitation file input source directory:  '//precdir
      end if
    else                         !just set results to something
      dtprec=0._wp
      precdir=''
    end if

    !ELECTRIC FIELD FILE INPUT INFORMATION
    read(u,*) flagE0file
    if (flagE0file==1) then    !get the location of the precipitation input files
      read(u,*) dtE0
    
      read(u,'(a256)') buf
      E0dir = trim(buf)
    
      if (myid==0) then
        print *, 'Electric field file input cadence (s):  ',dtE0
        print *, 'Electric field file input source directory:  '//E0dir
      end if
    else                         !just set results to something
      dtE0=0._wp
      E0dir=''
    end if

    close(u)
  end subroutine read_configfile


  subroutine create_outdir(outdir,infile,indatsize,indatgrid,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)

    !------------------------------------------------------------
    !-------CREATES OUTPUT DIRECTORY, MOVES CONFIG FILES THERE AND
    !-------GENERATES A GRID OUTPUT FILE
    !------------------------------------------------------------
    
    character(*), intent(in) :: outdir, & !command line argument output directory
                                infile, & !command line argument input file
                                indatsize,indatgrid,sourcedir, precdir,E0dir
    integer, intent(in) :: flagdneu, flagprecfile, flagE0file
    
    integer :: ierr

    !MAKE A COPY OF THE INPUT DATA IN THE OUTPUT DIRECTORY (MAYBE SHOULD COPY SOURCE CODE TOO???)
    call execute_command_line('mkdir -pv '//outdir//'/inputs', exitstat=ierr)
    if (ierr /= 0) error stop 'error creating output directory' 

    call execute_command_line('cp -rv '//infile//' '//outdir//'/inputs/', exitstat=ierr)
    if (ierr /= 0) error stop 'error copying input parameters to output directory' 
    call execute_command_line('cp -rv '//indatsize//' '//outdir//'/inputs/', exitstat=ierr)
    if (ierr /= 0) error stop 'error copying input parameters to output directory' 
    call execute_command_line('cp -rv '//indatgrid//' '//outdir//'/inputs/', exitstat=ierr)
    if (ierr /= 0) error stop 'error copying input parameters to output directory' 
    call execute_command_line('cp -rv '//indatfile//' '//outdir//'/inputs/', exitstat=ierr)
    if (ierr /= 0) error stop 'error copying input parameters to output directory' 

    !MAKE COPIES OF THE INPUT DATA, AS APPROPRIATE
    if (flagdneu/=0) then
      call execute_command_line('mkdir -pv '//outdir//'/inputs/neutral_inputs')
      call execute_command_line('cp -rv '//sourcedir//'/* '//outdir//'/inputs/neutral_inputs/', exitstat=ierr)
    end if
    if (ierr /= 0) error stop 'error copying neutral input parameters to output directory' 

    if (flagprecfile/=0) then
      call execute_command_line('mkdir -pv '//outdir//'/inputs/prec_inputs')
      call execute_command_line('cp -rv '//precdir//'/* '//outdir//'/inputs/prec_inputs/', exitstat=ierr)
    end if
    if (ierr /= 0) error stop 'error copying input precipitation parameters to output directory' 
    
    if (flagE0file/=0) then
      call execute_command_line('mkdir -pv '//outdir//'/inputs/Efield_inputs')
      call execute_command_line('cp -rv '//E0dir//'/* '//outdir//'/inputs/Efield_inputs/', exitstat=ierr)
    end if
    if (ierr /= 0) error stop 'error copying input energy parameters to output directory' 

    !NOW STORE THE VERSIONS/COMMIT IDENTIFIER IN A FILE IN THE OUTPUT DIRECTORY
    ! this can break on POSIX due to copying files in endless loop, commended out - MH
    !call execute_command_line('mkdir -pv '//outdir//'/inputs/source/', exitstat=ierr)
    !if (ierr /= 0) error stop 'error creating input source parameter output directory'
    !call execute_command_line('cp -rv ./* '//outdir//'/inputs/source/', exitstat=ierr)
    !if (ierr /= 0) error stop 'error creating input source parameter output directory' 
    
    call execute_command_line('git rev-parse --short HEAD > '//outdir//'/gitrev.log')

  end subroutine create_outdir


  subroutine create_outdir_mag(outdir,fieldpointfile)

    !------------------------------------------------------------
    !-------CREATES OUTPUT DIRECTORY FOR MAGNETIC FIELD CALCULATIONS
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    character(*), intent(in) :: fieldpointfile


    !NOTE HERE THAT WE INTERPRET OUTDIR AS THE BASE DIRECTORY CONTAINING SIMULATION OUTPUT
    call execute_command_line('mkdir '//outdir//'/magfields/')
    call execute_command_line('mkdir '//outdir//'/magfields/input/')
    call execute_command_line('cp '//fieldpointfile//' '//outdir//'/magfields/input/magfieldpoints.dat')

  end subroutine create_outdir_mag


  subroutine input_plasma(x1,x2,x3all,indatsize,ns,vs1,Ts)

    !------------------------------------------------------------
    !-------A BASIC WRAPPER FOR THE ROOT AND WORKER INPUT FUNCTIONS
    !-------BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
    !-------VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)
    !------------------------------------------------------------

    real(8), dimension(-1:), intent(in) :: x1, x2, x3all
    character(*), intent(in) :: indatsize

    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts

 
    if (myid==0) then
      !ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
      print *, 'Assembling initial condition on root using '//indatsize//' '//indatfile
      call input_root_mpi(x1,x2,x3all,indatsize,ns,vs1,Ts)
    else
      !WORKERS RECEIVE THE IC DATA FROM ROOT
      call input_workers_mpi(ns,vs1,Ts)
    end if

  end subroutine input_plasma


  subroutine input_workers_mpi(ns,vs1,Ts)
 
    !------------------------------------------------------------
    !-------RECEIVE INITIAL CONDITIONS FROM ROOT PROCESS
    !------------------------------------------------------------
 
    real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts     
  
    call bcast_recv(ns,tagns)
    call bcast_recv(vs1,tagvs1)
    call bcast_recv(Ts,tagTs)
  
  end subroutine input_workers_mpi
  
  
  subroutine input_root_mpi(x1,x2,x3all,indatsize,ns,vs1,Ts)
 
    !------------------------------------------------------------
    !-------READ INPUT FROM FILE AND DISTRIBUTE TO WORKERS.  
    !-------STATE VARS ARE EXPECTED INCLUDE GHOST CELLS.  NOTE ALSO
    !-------THAT RECORD-BASED INPUT IS USED SO NO FILES > 2GB DUE
    !-------TO GFORTRAN BUG WHICH DISALLOWS 8 BYTE INTEGER RECORD
    !-------LENGTHS.
    !------------------------------------------------------------
    
    real(8), dimension(-1:), intent(in) :: x1, x2, x3all
    character(*), intent(in) :: indatsize
    real(8), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts    

    integer :: lx1,lx2,lx3,lx3all,isp
    
    real(8), dimension(-1:size(x1,1)-2,-1:size(x2,1)-2,-1:size(x3all,1)-2,1:lsp) :: nsall,vs1all,Tsall
    real(8), dimension(:,:,:,:), allocatable :: statetmp
    integer :: lx1in,lx2in,lx3in,u
    real(8) :: tin
    real(8), dimension(3) :: ymdtmp
    
    real(8) :: tstart,tfin    

    !SYSTEM SIZES
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4
    lx3all=size(nsall,3)-4
    
            
    !READ IN FROM FILE, AS OF CURVILINEAR BRANCH THIS IS NOW THE ONLY INPUT
    !OPTION
    open(newunit=u,file=indatsize,status='old',form='unformatted', access='stream', action='read')
    read(u) lx1in,lx2in,lx3in
    close(u)
    print *, 'Input file has size:  ',lx1in,lx2in,lx3in
    print *, 'Target grid structure has size',lx1,lx2,lx3all

    if (flagswap==1) then
      print *, '2D simulations grid detected, swapping input file dimension sizes and permuting input arrays'
      lx3in=lx2in
      lx2in=1
    end if

    if (.not. (lx1==lx1in .and. lx2==lx2in .and. lx3all==lx3in)) then
      error stop '!!!The input data must be the same size as the grid which you are running the simulation on' // & 
           '- use a script to interpolate up/down to the simulation grid'
    end if

    open(newunit=u,file=indatfile,status='old',form='unformatted', access='stream', action='read')
    read(u) ymdtmp,tin

    if (flagswap/=1) then
      read(u) nsall(1:lx1,1:lx2,1:lx3all,1:lsp)
      read(u) vs1all(1:lx1,1:lx2,1:lx3all,1:lsp)
      read(u) Tsall(1:lx1,1:lx2,1:lx3all,1:lsp)
    else
      allocate(statetmp(lx1,lx3all,lx2,lsp))
      read(u) statetmp
      
      print *, shape(statetmp),shape(nsall)
      nsall(1:lx1,1:lx2,1:lx3all,1:lsp)=reshape(statetmp,[lx1,lx2,lx3all,lsp],order=[1,3,2,4])
      read(u) statetmp
      vs1all(1:lx1,1:lx2,1:lx3all,1:lsp)=reshape(statetmp,[lx1,lx2,lx3all,lsp],order=[1,3,2,4])

      read(u) statetmp
      Tsall(1:lx1,1:lx2,1:lx3all,1:lsp)=reshape(statetmp,[lx1,lx2,lx3all,lsp],order=[1,3,2,4])    !permute the dimensions so that 2D runs are parallelized
      deallocate(statetmp)
    end if
    close(u)
    print *, 'Done gathering input...'


    !USER SUPPLIED FUNCTION TO TAKE A REFERENCE PROFILE AND CREATE INITIAL CONDITIONS FOR ENTIRE GRID.  ASSUMING THAT THE
    !INPUT DATA ARE EXACTLY THE CORRECT SIZE (AS IS THE CASE WITH FILE INPUT) THIS IS NOW SUPERFLUOUS
    print *, 'Done setting initial conditions...'
    print *, 'Min/max input density:  ',minval(pack(nsall(:,:,:,7),.true.)),maxval(pack(nsall(:,:,:,7),.true.))
    print *, 'Min/max input velocity:  ',minval(pack(vs1all(:,:,:,:),.true.)),maxval(pack(vs1all(:,:,:,:),.true.))
    print *, 'Min/max input temperature:  ',minval(pack(Tsall(:,:,:,:),.true.)),maxval(pack(Tsall(:,:,:,:),.true.))


    !ROOT BROADCASTS IC DATA TO WORKERS
    call cpu_time(tstart)
    call bcast_send(nsall,tagns,ns)
    call bcast_send(vs1all,tagvs1,vs1)
    call bcast_send(Tsall,tagTs,Ts)
    call cpu_time(tfin)
    print *, 'Done sending ICs to workers...  Time:  ',tfin-tstart

  end subroutine input_root_mpi


  subroutine input_plasma_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)

    !------------------------------------------------------------
    !-------READS, AS INPUT, A FILE GENERATED BY THE GEMINI.F90 PROGRAM.
    !-------THIS SUBROUTINE IS A WRAPPER FOR SEPARATE ROOT/WORKER
    !-------CALLS
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec
    real(8), dimension(:,:,:), intent(out) :: J1,J2,J3


    if (myid==0) then
      !ROOT FINDS/CALCULATES INITIAL CONDITIONS AND SENDS TO WORKERS
      print *, 'Assembling current density data on root...  '
      call input_root_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
    else
      !WORKERS RECEIVE THE IC DATA FROM ROOT
      call input_workers_currents(J1,J2,J3)
    end if

  end subroutine input_plasma_currents


  subroutine input_root_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)

    !------------------------------------------------------------
    !-------READS, AS INPUT, A FILE GENERATED BY THE GEMINI.F90 PROGRAM
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput
    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec
    real(8), dimension(:,:,:), intent(out) :: J1,J2,J3

    real(8), dimension(:,:,:), allocatable :: tmparray3D
    real(8), dimension(:,:,:,:), allocatable :: tmparray4D
    character(:), allocatable :: filenamefull
    real(8), dimension(:,:,:), allocatable :: J1all,J2all,J3all
    real(8), dimension(:,:,:), allocatable :: tmpswap
    real(8) :: tmpdate
    integer :: u


    !CHECK TO MAKE SURE WE ACTUALLY HAVE THE DATA WE NEED TO DO THE MAG COMPUTATIONS.
    if (flagoutput==3) error stop '  !!!I need current densities in the output to compute magnetic fields!!!'


    !FORM THE INPUT FILE NAME
    filenamefull=date_filename(outdir,ymd,UTsec)
    print *, 'Input file name for current densities:  ',filenamefull
    open(newunit=u,file=filenamefull,status='old',form='unformatted',access='stream',action='read')
    read(u) tmpdate
    write(*,*) 'File year:  ',tmpdate
    read(u) tmpdate
    write(*,*) 'File month:  ',tmpdate
    read(u) tmpdate
    write(*,*) 'File day:  ',tmpdate
    read(u) tmpdate
    write(*,*) 'File UThrs:  ',tmpdate


    !LOAD THE DATA
    if (flagoutput==2) then    !the simulation data have only averaged plasma parameters
      write(*,*) '  Reading in files containing averaged plasma parameters of size:  ',lx1*lx2*lx3all
      allocate(tmparray3D(lx1,lx2,lx3all))
      !MZ:  I've found what I'd consider to be a gfortran bug here.  If I read
      !in a flat array (i.e. a 1D set of data) I hit EOF, according to runtime
      !error, well before I'm actually out of data this happens with a 20GB
      !input file for not for a 3GB input file...  This doesn't happen when I do
      !the reading with 3D arrays.  
      read(u) tmparray3D    !ne - could be done with some judicious fseeking...
      read(u) tmparray3D    !vi
      read(u) tmparray3D    !Ti
      read(u) tmparray3D    !Te
      deallocate(tmparray3D)
    else    !full output parameters are in the output files
      write(*,*) '  Reading in files containing full plasma parameters of size:  ',lx1*lx2*lx3all*lsp
      allocate(tmparray4D(lx1,lx2,lx3all,lsp))
      read(u) tmparray4D
      read(u) tmparray4D
      read(u) tmparray4D
      read(u) tmparray4D
      deallocate(tmparray4D)
    end if


    !PERMUTE THE ARRAYS IF NECESSARY
    write(*,*) '  File fast-forward done, now reading currents...'
    allocate(J1all(lx1,lx2,lx3all),J2all(lx1,lx2,lx3all),J3all(lx1,lx2,lx3all))
    if (flagswap==1) then
      allocate(tmpswap(lx1,lx3all,lx2))
      read(u) tmpswap
      J1all=reshape(tmpswap,[lx1,lx2,lx3all],order=[1,3,2])
      read(u) tmpswap
      J2all=reshape(tmpswap,[lx1,lx2,lx3all],order=[1,3,2])
      read(u) tmpswap
      J3all=reshape(tmpswap,[lx1,lx2,lx3all],order=[1,3,2])
      deallocate(tmpswap)
    else    !no need to permute dimensions for 3D simulations
      read(u) J1all,J2all,J3all
    end if
    write(*,*) 'Min/max current data:  ',minval(J1all),maxval(J1all),minval(J2all),maxval(J2all),minval(J3all),maxval(J3all)


    !DISTRIBUTE DATA TO WORKERS AND TAKE A PIECE FOR ROOT
    call bcast_send(J1all,tagJ1,J1)
    call bcast_send(J2all,tagJ2,J2)
    call bcast_send(J3all,tagJ3,J3)


    !CLEAN UP MEMORY
    deallocate(J1all,J2all,J3all)

  end subroutine input_root_currents


  subroutine input_workers_currents(J1,J2,J3)

    !------------------------------------------------------------
    !-------WORKER INPUT FUNCTIONS FOR GETTING CURRENT DENSITIES
    !------------------------------------------------------------

    real(8), dimension(:,:,:), intent(out) :: J1,J2,J3


    !ALL WE HAVE TO DO IS WAIT TO RECEIVE OUR PIECE OF DATA FROM ROOT
    call bcast_recv(J1,tagJ1)
    call bcast_recv(J2,tagJ2)
    call bcast_recv(J3,tagJ3)

  end subroutine input_workers_currents


  subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)

    !------------------------------------------------------------
    !-------A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
    !-------BOTH ROOT AND WORKERS CALL THIS PROCEDURE SO UNALLOCATED
    !-------VARIABLES MUST BE DECLARED AS ALLOCATABLE, INTENT(INOUT)
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput

    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts

    real(8), dimension(:,:,:), allocatable, intent(inout) :: Phiall     !these jokers may not be allocated, but this is allowed as of f2003
    real(8), dimension(:,:,:), intent(in) :: J1,J2,J3


    if (myid/=0) then
      call output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3)
    else
      call output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)    !only option that works with >2GB files
    end if  

  end subroutine output_plasma


  subroutine output_workers_mpi(vs2,vs3,ns,vs1,Ts,J1,J2,J3)

    !------------------------------------------------------------
    !-------SEND COMPLETE DATA FROM WORKERS TO ROOT PROCESS FOR OUTPUT.  
    !-------STATE VARS ARE EXPECTED TO INCLUDE GHOST CELLS
    !------------------------------------------------------------
  
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts     
    real(8), dimension(:,:,:), intent(in) :: J1,J2,J3

    integer :: lx1,lx2,lx3,lx3all,isp
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg


    !SYSTEM SIZES (W/O GHOST CELLS)
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4
    
    
    !ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
    v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
    v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)
  

    !SEND MY GRID DATA TO THE ROOT PROCESS
    call gather_send(v2avg,tagv2)
    call gather_send(v3avg,tagv3)
    call gather_send(ns,tagns)
    call gather_send(vs1,tagvs1)
    call gather_send(Ts,tagTs)


    !------- SEND ELECTRODYNAMIC PARAMETERS TO ROOT
    call gather_send(J1,tagJ1)
    call gather_send(J2,tagJ2)
    call gather_send(J3,tagJ3)  

  end subroutine output_workers_mpi


  subroutine output_root_stream_mpi(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)

    !------------------------------------------------------------
    !-------COLLECT OUTPUT FROM WORKERS AND WRITE TO A FILE USING
    !-------STREAM I/O.    
    !-------STATE VARS ARE EXPECTED INCLUDE GHOST CELLS
    !------------------------------------------------------------
  
    character(*), intent(in) :: outdir
    integer, intent(in) :: flagoutput

    integer, dimension(3), intent(in) :: ymd
    real(8), intent(in) :: UTsec
    real(8), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts    
    
    real(8), dimension(:,:,:), intent(in) :: Phiall
    real(8), dimension(:,:,:), intent(in) :: J1,J2,J3
 
    integer :: lx1,lx2,lx3,lx3all,isp, u
    real(8), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: v2avg,v3avg
    real(8), dimension(-1:size(Phiall,1)+2,-1:size(Phiall,2)+2,-1:size(Phiall,3)+2,1:lsp) :: nsall,vs1all,Tsall
    real(8), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: v2avgall,v3avgall,v1avgall,Tavgall,neall,Teall
    real(8), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: J1all,J2all,J3all
    character(:), allocatable :: filenamefull
    integer(8) :: recordlength   !can be 8 byte with compiler flag -frecord-marker=8

    real(8), dimension(:,:,:), allocatable :: permarray,tmparray    !permuted variables to be allocated for 2D output  

 
    !SYSTEM SIZES
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4
    lx3all=size(Phiall,3)
    
    
    !ONLY AVERAGE DRIFTS PERP TO B NEEDED FOR OUTPUT
    v2avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs2(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v2avg=v2avg/ns(1:lx1,1:lx2,1:lx3,lsp)    !compute averages for output.
    v3avg=sum(ns(1:lx1,1:lx2,1:lx3,1:lsp-1)*vs3(1:lx1,1:lx2,1:lx3,1:lsp-1),4)
    v3avg=v3avg/ns(1:lx1,1:lx2,1:lx3,lsp)


    !GET THE SUBGRID DATA FORM THE WORKERS  
    call gather_recv(v2avg,tagv2,v2avgall) 
    call gather_recv(v3avg,tagv3,v3avgall)
    call gather_recv(ns,tagns,nsall)
    call gather_recv(vs1,tagvs1,vs1all)
    call gather_recv(Ts,tagTs,Tsall)


    !RADD--- NEED TO ALSO GATHER FULL GRID ELECTRODYANMICS PARAMTERS FROM WORKERS
    call gather_recv(J1,tagJ1,J1all)
    call gather_recv(J2,tagJ2,J2all)
    call gather_recv(J3,tagJ3,J3all)


    !COMPUTE AVERAGE VALUES FOR ION PLASMA PARAMETERS
    v1avgall=sum(nsall(1:lx1,1:lx2,1:lx3all,1:lsp-1)*vs1all(1:lx1,1:lx2,1:lx3all,1:lsp-1),4)
    v1avgall=v1avgall/nsall(1:lx1,1:lx2,1:lx3all,lsp)    !compute averages for output.
    Tavgall=sum(nsall(1:lx1,1:lx2,1:lx3all,1:lsp-1)*Tsall(1:lx1,1:lx2,1:lx3all,1:lsp-1),4)
    Tavgall=Tavgall/nsall(1:lx1,1:lx2,1:lx3all,lsp)    !compute averages for output.
    neall=nsall(1:lx1,1:lx2,1:lx3all,lsp)
    Teall=Tsall(1:lx1,1:lx2,1:lx3all,lsp)


    !FIGURE OUT THE FILENAME
    filenamefull=date_filename(outdir,ymd,UTsec)
    print *, 'Output file name:  ',filenamefull


    !SOME DEBUG OUTPUT ON FILE SIZE
    recordlength=int(8,8)+int(8,8)*int(3,8)*int(lx1,8)*int(lx2,8)*int(lx3all,8)*int(lsp,8)+ &
                 int(8,8)*int(5,8)*int(lx1,8)*int(lx2,8)*int(lx3all,8)+ &
                 int(8,8)*int(lx2,8)*int(lx3all,8)
    print *, 'Output bit length:  ',recordlength,lx1,lx2,lx3all,lsp


    !WRITE THE DATA
    open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')    !has no problem with > 2GB output files
    write(u) real(ymd,wp),UTsec/3600._wp    !no matter what we must output date and time

    if (flagswap/=1) then
      select case (flagoutput)
        case (2)    !output ISR-like average parameters
          write(u) neall(1:lx1,1:lx2,1:lx3all),v1avgall(1:lx1,1:lx2,1:lx3all), &    !output of ISR-like parameters (ne,Ti,Te,v1,etc.)
                      Tavgall(1:lx1,1:lx2,1:lx3all),Teall(1:lx1,1:lx2,1:lx3all),J1all(1:lx1,1:lx2,1:lx3all), &
                      J2all(1:lx1,1:lx2,1:lx3all), &
                      J3all(1:lx1,1:lx2,1:lx3all),v2avgall(1:lx1,1:lx2,1:lx3all),v3avgall(1:lx1,1:lx2,1:lx3all)
        case (3)     !just electron density
          print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
          write(u) neall(1:lx1,1:lx2,1:lx3all)
        case default    !output everything
          print *, '!!!NOTE:  Input file has selected full ouptut, large files may result!'
          write(u) nsall(1:lx1,1:lx2,1:lx3all,:),vs1all(1:lx1,1:lx2,1:lx3all,:), &    !this is full output of all parameters in 3D
                      Tsall(1:lx1,1:lx2,1:lx3all,:),J1all(1:lx1,1:lx2,1:lx3all),J2all(1:lx1,1:lx2,1:lx3all), &
                      J3all(1:lx1,1:lx2,1:lx3all),v2avgall(1:lx1,1:lx2,1:lx3all),v3avgall(1:lx1,1:lx2,1:lx3all)
        end select
    else                 !2D simulation for which arrays were permuted
      print *, '!!!NOTE:  Permuting arrays prior to output...'
      select case (flagoutput)
        case (2)    !averaged parameters
          allocate(permarray(lx1,lx3all,lx2))    !temporary work array that has been permuted
          permarray=reshape(neall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(v1avgall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(Tavgall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(Teall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(J1all,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(J3all,[lx1,lx3all,lx2],order=[1,3,2])    !Note that components need to be swapped too
          write(u) permarray
          permarray=reshape(J2all,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(v3avgall,[lx1,lx3all,lx2],order=[1,3,2])    !Note swapping of components
          write(u) permarray
          permarray=reshape(v2avgall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          deallocate(permarray) 
        case (3)     !electron density only output
          print *, '!!!NOTE:  Input file has selected electron density only output, make sure this is what you really want!'
          allocate(permarray(lx1,lx3all,lx2))    !temporary work array that has been permuted
          permarray=reshape(neall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          deallocate(permarray)
        case default
          print *, '!!!NOTE:  Input file has selected full ouptut, large files may result!'
          allocate(permarray(lx1,lx3all,lx2))    !temporary work array that has been permuted
          allocate(tmparray(lx1,lx2,lx3all))
          do isp=1,lsp
            tmparray=nsall(1:lx1,1:lx2,1:lx3all,isp)
            permarray=reshape(tmparray,[lx1,lx3all,lx2],order=[1,3,2])
            write(u) permarray
          end do
          do isp=1,lsp
            tmparray=vs1all(1:lx1,1:lx2,1:lx3all,isp)
            permarray=reshape(tmparray,[lx1,lx3all,lx2],order=[1,3,2])
            write(u) permarray
          end do 
          do isp=1,lsp
            tmparray=Tsall(1:lx1,1:lx2,1:lx3all,isp)
            permarray=reshape(tmparray,[lx1,lx3all,lx2],order=[1,3,2])
            write(u) permarray
          end do
          permarray=reshape(J1all,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(J3all,[lx1,lx3all,lx2],order=[1,3,2])    !Note that components need to be swapped too
          write(u) permarray
          permarray=reshape(J2all,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          permarray=reshape(v3avgall,[lx1,lx3all,lx2],order=[1,3,2])    !Note swapping of components
          write(u) permarray
          permarray=reshape(v2avgall,[lx1,lx3all,lx2],order=[1,3,2])
          write(u) permarray
          deallocate(permarray)
          deallocate(tmparray)
      end select 
    end if
    if (gridflag==1) then
      print *, 'Writing topside boundary conditions for inverted-type grid...'
      write(u)  Phiall(1,:,:)
    else
      print *, 'Writing topside boundary conditions for non-inverted-type grid...'
      write(u)  Phiall(lx1,:,:)
    end if

    close(u)
  
  end subroutine output_root_stream_mpi


  subroutine output_magfields(outdir,ymd,UTsec,Br,Btheta,Bphi)

    !------------------------------------------------------------
    !-------A BASIC WRAPPER FOR THE ROOT AND WORKER OUTPUT FUNCTIONS
    !-------BOTH ROOT AND WORKERS CALL THIS PROCEDURE TO GENERATE
    !-------MAGNETIC FIELD OUTPUT FILES,  WE ASSUME THE ROOT PROCESS
    !-------HAS ALREADY REDUCED THE MAGNETIC FIELD DATA
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, intent(in) :: ymd(3)
    real(8), intent(in) :: UTsec
    real(8), dimension(:), intent(in)  :: Br,Btheta,Bphi

    character(:), allocatable :: outdir_composite, filenamefull
    integer :: u


    !FORM THE INPUT FILE NAME
    outdir_composite=outdir//'/magfields/'
    filenamefull=date_filename(outdir_composite,ymd,UTsec)
    print *, '  Output file name (magnetic fields):  ',filenamefull
    open(newunit=u,file=filenamefull,status='replace',form='unformatted',access='stream',action='write')


    !DUMP THE OUTPUT DATA
    write(u) Br,Btheta,Bphi 
    close(u)
  end subroutine output_magfields


  pure function date_filename(outdir,ymd,UTsec)

    !------------------------------------------------------------
    !-------GENERATE A FILENAME STRING OUT OF A GIVEN DATE/TIME
    !------------------------------------------------------------

    character(*), intent(in) :: outdir
    integer, intent(in) :: ymd(3)
    real(8), intent(in) :: UTsec
    character(:), allocatable :: date_filename

    character(16) :: ssec
    character(9) :: symd


    ! UTC second (float, 0.0 .. 86400) 
    write(ssec,'(f12.6,a4)') UTsec,'.dat'    !file name that has 6 decimal points on time stamp

    ! year_month_day
    write(symd,'(i4,2i0.2,a1)') ymd, '_'
    
    ! assemble
    date_filename = outdir // '/' // symd // ssec

  end function date_filename

end module io
