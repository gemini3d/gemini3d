module potential_comm

!------------------------------------------------------------
!-------NOTE: THIS MODULE, WHICH HANDLES COMMUNICATION FOR
!-------MESSAGE PASSING PREP OF POTENTIAL SOLVER IS ISOLATED
!-------FROM THE FUNCTION WHICH CALLS ILUPACK DUE TO THE
!-------APPARENT NEED FOR ILU CALLING FUNCTIONS TO DEFAULT TO
!-------8-BYTE INTEGER TYPES.  THIS WILL APPARENTLY NOT WORK FOR ANY
!-------MPI PROGRAM, WHICH USES 4 BYTE INTEGERS COMMON BLOCKS.
!-------
!-------FOR NOW THIS MODULE IGNORES POLARZATION CURRENT EFFECTS.
!------------------------------------------------------------


use phys_consts, only : pi
use grid, only : curvmesh
use calculus
use collisions
use potential_ilupack, only : elliptic3D,elliptic3D_pol_conv
use mpimod
implicit none


!INTERFACE REQUIRED FOR ASSUMED-SHAPE ARGUMENTS
interface
  subroutine potentialBCs3D(t,sig0all,x1,x2,x3all,dx2,dx3all,Vminx1,Vmaxx1, &
                      Vminx2,Vmaxx2,Vminx3,Vmaxx3,E01all,E02all,E03all,flagdirich)

    real(wp), intent(in) :: t
    real(wp), dimension(:,:,:), intent(in) ::  sig0all
    real(wp), dimension(-1:), intent(in) :: x1
    real(wp), dimension(-1:), intent(in) :: x2
    real(wp), dimension(-1:), intent(in) :: x3all
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(0:), intent(in) :: dx3all

    real(wp), dimension(:,:), intent(out) :: Vminx1,Vmaxx1
    real(wp), dimension(:,:), intent(out) :: Vminx2,Vmaxx2
    real(wp), dimension(:,:), intent(out) :: Vminx3,Vmaxx3
    real(wp), dimension(:,:,:), intent(out) :: E01all,E02all,E03all
    integer(8), intent(out) :: flagdirich

  end subroutine potentialBCs3D
end interface


!OVERLOAD THE SOLVERS FOR CURVILINEAR MESHES
interface electrodynamics
 module procedure electrodynamics_cart,electrodynamics_curv
end interface electrodynamics

contains


  subroutine electrodynamics_cart(it,t,dt,nn,Tn,ns,Ts,vs1,B1,vs2,vs3, &
                             x1,x2,x3all,dx1,dx1i,dx2,dx2i,dx3,dx3i,dx3all,dx3iall, &
                             potsolve, &
                             E1,E2,E3,J1,J2,J3, &
                             E1all,E2all,E3all,J1all,J2all,J3all,Phiall)

    !------------------------------------------------------------
    !-------THIS IS A WRAPPER FUNCTION FOR THE ELECTRODYANMICS
    !-------PART OF THE MODEL.  BOTH THE ROOT AND WORKER PROCESSES
    !-------CALL THIS SAME SUBROUTINE, WHEN THEN BRANCHES INTO
    !-------DIFFERENT TASKS FOR EACH AFTER ALL COMPUTE CONDUCTIVITIES
    !-------AND INERTIAL CAPACITANCE.
    !-------
    !-------NOTE THAT THE ALLOCATION STATUS
    !-------OF THE ALL VARIABLES FOR THE WORKERS WILL BE UN-
    !-------ALLOCATED.  SO THIS CODE IS ONLY COMPLIANT WITH
    !-------FORTRAN >= 2003 STANDARDS.
    !------------------------------------------------------------

    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt

    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1

    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3

    real(wp), dimension(-1:), intent(in) :: x1
    real(wp), dimension(-1:), intent(in) :: x2
    real(wp), dimension(-1:), intent(in) :: x3all   !workers actually are assumed to have a copy of grid all variable, otherwise and allocatable is needed here
    real(wp), dimension(0:), intent(in) :: dx1
    real(wp), dimension(:), intent(in) :: dx1i
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(:), intent(in) :: dx2i
    real(wp), dimension(0:), intent(in) :: dx3
    real(wp), dimension(:), intent(in) :: dx3i
    real(wp), dimension(0:), intent(in) :: dx3all
    real(wp), dimension(:), intent(in) :: dx3iall

    integer, intent(in) :: potsolve

    real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: E1all,E2all,E3all
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: J1all,J2all,J3all
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall

    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: sig0,sigP,sigH
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp) :: muP,muH
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: incap
    real(wp) :: tstart,tfin

    integer :: lx1,lx2,lx3,isp


    !SIZES
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4


    !POTENTIAL SOLUTION (IF REQUIRED)
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH)
    call capacitance(ns,B1,incap)
  !  incap=0d0   !will get rid of polarization terms without having to switch solver
    if (potsolve == 1 .or. potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
      call cpu_time(tstart)

      if (myid/=0) then    !role-specific communication pattern (all-to-root-to-all), workers initiate with sends
        call potential_workers_mpi(sig0,sigP,sigH,incap,vs2,vs3,potsolve, &
                                   E1,E2,E3,J1,J2,J3)
      else
        call potential_root_mpi_cart(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3, &
                                  x1,x2,x3all,dx1,dx1i,dx2,dx2i,dx3,dx3i,dx3all,dx3iall, &
                                  potsolve, &
                                  E1,E2,E3,J1,J2,J3, &
                                  E1all,E2all,E3all,J1all,J2all,J3all,Phiall)
      end if

      do isp=1,lsp
        vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3
        vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3
      end do

  !    do isp=1,lsp    !To leading order the ion drifts do not include the polarization parts, otherwise it may mess up polarization convective term in the electrodynamics solver...
  !      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
  !      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
  !    end do

      call cpu_time(tfin)

      if (myid==0) then
        write(*,*) 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
      end if
    else if (potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
      !Do nothing for now...
    else   !null solve
      E1=0d0; E2=0d0; E3=0d0; J1=0d0; J2=0d0; J3=0d0;
      vs2=0d0; vs3=0d0;
    end if

  end subroutine electrodynamics_cart


  subroutine electrodynamics_curv(it,t,dt,nn,Tn,ns,Ts,vs1,B1,vs2,vs3,x, &
                             potsolve, &
                             E1,E2,E3,J1,J2,J3, &
                             E1all,E2all,E3all,J1all,J2all,J3all,Phiall)

    !------------------------------------------------------------
    !-------THIS IS A WRAPPER FUNCTION FOR THE ELECTRODYANMICS
    !-------PART OF THE MODEL.  BOTH THE ROOT AND WORKER PROCESSES
    !-------CALL THIS SAME SUBROUTINE, WHEN THEN BRANCHES INTO
    !-------DIFFERENT TASKS FOR EACH AFTER ALL COMPUTE CONDUCTIVITIES
    !-------AND INERTIAL CAPACITANCE.
    !-------
    !-------NOTE THAT THE ALLOCATION STATUS
    !-------OF THE ALL VARIABLES FOR THE WORKERS WILL BE UN-
    !-------ALLOCATED.  SO THIS CODE IS ONLY COMPLIANT WITH
    !-------FORTRAN >= 2003 STANDARDS.
    !------------------------------------------------------------

    integer, intent(in) :: it
    real(wp), intent(in) :: t,dt

    real(wp), dimension(:,:,:,:), intent(in) :: nn
    real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts,vs1
    real(wp), dimension(-1:,-1:,-1:), intent(in) :: B1

    real(wp), dimension(-1:,-1:,-1:,:), intent(inout) ::  vs2,vs3

    type(curvmesh), intent(in) :: x

    integer, intent(in) :: potsolve

    real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: E1all,E2all,E3all
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: J1all,J2all,J3all
    real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall

    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: sig0,sigP,sigH
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4,lsp) :: muP,muH
    real(wp), dimension(1:size(ns,1)-4,1:size(ns,2)-4,1:size(ns,3)-4) :: incap
    real(wp) :: tstart,tfin

    integer :: lx1,lx2,lx3,isp


    !SIZES
    lx1=size(ns,1)-4
    lx2=size(ns,2)-4
    lx3=size(ns,3)-4


    !POTENTIAL SOLUTION (IF REQUIRED)
    call conductivities(nn,Tn,ns,Ts,vs1,B1,sig0,sigP,sigH,muP,muH)
    call capacitance(ns,B1,incap)
  !  incap=0d0   !will get rid of polarization terms without having to switch solver
    if (potsolve == 1 .or. potsolve == 3) then    !electrostatic solve or electrostatic alt. solve
      call cpu_time(tstart)

      if (myid/=0) then    !role-specific communication pattern (all-to-root-to-all), workers initiate with sends
        call potential_workers_mpi(sig0,sigP,sigH,incap,vs2,vs3,potsolve, &
                                   E1,E2,E3,J1,J2,J3)
      else
        call potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,x, &
                                  potsolve, &
                                  E1,E2,E3,J1,J2,J3, &
                                  E1all,E2all,E3all,J1all,J2all,J3all,Phiall)
      end if

      do isp=1,lsp
        vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3
        vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3
      end do

  !    do isp=1,lsp    !To leading order the ion drifts do not include the polarization parts, otherwise it may mess up polarization convective term in the electrodynamics solver...
  !      vs2(1:lx1,1:lx2,1:lx3,isp)=muP(:,:,:,isp)*E2-muH(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE2Dt
  !      vs3(1:lx1,1:lx2,1:lx3,isp)=muH(:,:,:,isp)*E2+muP(:,:,:,isp)*E3+ms(isp)/qs(isp)/B1**2*DE3Dt
  !    end do

      call cpu_time(tfin)

      if (myid==0) then
        write(*,*) 'Potential solution for time step:  ',t,' took ',tfin-tstart,' seconds...'
      end if
    else if (potsolve == 2) then  !inductive form of model, could this be subcycled to speed things up?
      !Do nothing for now...
    else   !null solve
      E1=0d0; E2=0d0; E3=0d0; J1=0d0; J2=0d0; J3=0d0;
      vs2=0d0; vs3=0d0;
    end if

  end subroutine electrodynamics_curv


  subroutine potential_root_mpi_cart(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3, &
                                x1,x2,x3all,dx1,dx1i, &
                                dx2,dx2i,dx3,dx3i,dx3all,dx3iall, &
                                potsolve, &
                                E1,E2,E3,J1,J2,J3, &
                                E1all,E2all,E3all,J1all,J2all,J3all,Phiall)

    !------------------------------------------------------------
    !-------ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL
    !------------------------------------------------------------

    integer, intent(in) :: it    !note that the makefile and ilupack insist that the default integer type is 8...
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH
    real(wp), dimension(:,:,:), intent(in) ::  incap
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3

    real(wp), dimension(-1:), intent(in) :: x1
    real(wp), dimension(-1:), intent(in) :: x2
    real(wp), dimension(-1:), intent(in) :: x3all
    real(wp), dimension(0:), intent(in) :: dx1
    real(wp), dimension(:), intent(in) :: dx1i
    real(wp), dimension(0:), intent(in) :: dx2
    real(wp), dimension(:), intent(in) :: dx2i
    real(wp), dimension(0:), intent(in) :: dx3
    real(wp), dimension(:), intent(in) :: dx3i
    real(wp), dimension(0:), intent(in) :: dx3all
    real(wp), dimension(:), intent(in) :: dx3iall

    integer, intent(in) :: potsolve

    real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
    real(wp), dimension(:,:,:), intent(inout) :: E1all,E2all,E3all
    real(wp), dimension(:,:,:), intent(out) :: J1all,J2all,J3all
    real(wp), dimension(:,:,:), intent(inout) :: Phiall   !not good form, but I'm lazy

    real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: sig0all,sigPall,sigHall   !really prefer sizes set by intent(in) vars.
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: incapall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: J1polall,J2polall,J3polall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: E1prevall,E2prevall,E3prevall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: DE2Dt,DE3Dt

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: srcterm
    real(wp), dimension(1:size(E1all,2),1:size(E1all,3)) :: Vminx1,Vmaxx1
    real(wp), dimension(1:size(E1all,1),1:size(E1all,3)) :: Vminx2,Vmaxx2
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2)) :: Vminx3,Vmaxx3
    integer(8) :: flagdirich    !ilupack expects 8-byte integers with the compile flags presently used

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: v2all,v3all   !stores drift velocs. for pol. current
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: E01all,E02all,E03all
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: DE2Dtall,DE3Dtall   !pol. drift

    integer :: ix2,ix3,lx1,lx2,lx3,lx3all

    real(wp) :: dVmaxx1,dE02all,dE03all

    integer :: iid,islstart,islfin


    !SIZES
    lx1=size(sig0,1)
    lx2=size(sig0,2)
    lx3=size(sig0,3)
    lx3all=size(E1all,3)

!
!    !COLLABORATE WITH WORKERS TO GET CONDUCTIVITIES
!    sigPall(:,:,1:lx3)=sigP(:,:,1:lx3)
!    sigHall(:,:,1:lx3)=sigH(:,:,1:lx3)
!    sig0all(:,:,1:lx3)=sig0(:,:,1:lx3)
!    incapall(:,:,1:lx3)=incap(:,:,1:lx3)
!    v2all(:,:,1:lx3)=vs2(1:lx1,1:lx2,1:lx3,1)    !use oxygen drift as ExB-ish drift...
!    v3all(:,:,1:lx3)=vs3(1:lx1,1:lx2,1:lx3,1)
!    do iid=1,lid-1
!      islstart=iid*lx3+1
!      islfin=islstart+lx3-1
!
!      call mpi_recv(sigPall(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagsigP,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!      call mpi_recv(sigHall(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagsigH,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!      call mpi_recv(sig0all(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagsig0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!      call mpi_recv(incapall(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagincap,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!      call mpi_recv(v2all(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagv2pol,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!      call mpi_recv(v3all(:,:,islstart:islfin),lx1*lx2*lx3,mpi_realprec, &
!                    iid,tagv3pol,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
!    end do
    !GATHER WORKER CONDUCTIVITY VALUES
    call gather_recv(sigP,tagsigP,sigPall)
    call gather_recv(sigH,tagsigH,sigHall)
    call gather_recv(sig0,tagsig0,sig0all)
    call gather_recv(incap,tagincap,incapall)
    v2=vs2(1:lx1,1:lx2,1:lx3,1)
    call gather_recv(v2,tagv2pol,v2all)
    v3=vs3(1:lx1,1:lx2,1:lx3,1)
    call gather_recv(v3,tagv3pol,v3all)


    !POPULATE BACKGROUND AND BOUNDARY CONDITION ARRAYS
    call potentialBCs3D(t,sig0all,x1,x2,x3all,dx2,dx3all,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                        E01all,E02all,E03all,flagdirich)
    J1all=0d0    !so this div is only perp components
    J2all=sigPall*E02all-sigHall*E03all    !BG x2 current
    J3all=sigHall*E02all+sigPall*E03all    !BG x3 current
    srcterm=div3D(J1all,J2all,J3all,dx1(1:lx1),dx2(1:lx2),dx3all(1:lx3all))


    !FORM AN INITIAL GUESS FOR THE ITERATIVE SOLVER AND CALL SOLVER
    if (it==1) then
      write(*,*) 'Using initial potential guess based on top BC...'
      if (flagdirich==1) then
        Phiall=spread(Vmaxx1,1,lx1)   !very good initial guess
      else
        Phiall=0d0     !will require many iterations to fix
      end if
    else
      write(*,*) 'Using initial potential guess based on previous time step solution...'
    end if


    !SOLVE FOR POTENTIAL USING INCOMPLETE LU FACTORIZATION
    if (mod(it,1)==0 .or. it==1) then    !Don't really need to evaluate every time, usually every 5 steps is okay...
      dVmaxx1=abs(maxval(pack(Vmaxx1,.true.))-minval(pack(Vmaxx1,.true.)))    !gfortran balks at massive if argument...
      dE02all=maxval(pack(abs(E02all),.true.))
      dE03all=maxval(pack(abs(E03all),.true.))

      if (dVmaxx1 < 1d-6 .and. dE02all < 1d-12 .and. dE03all < 1d-12) then   !also no eval if nothing is happening
        write(*,*) 'Potential boundary conditions AND source terms appear null, &
                    moving on with zero potential everywhere...'
        Phiall=0d0
      else
        write(*,*) 'Executing full solve for potential this time step...'
        Phiall=elliptic3D(srcterm,sig0all,sigPall,sigHall,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                        dx1,dx1i,dx2,dx2i,dx3all,dx3iall,Phiall,flagdirich)    !pure electrostatics, note all args are doubles,logicals so no integer type mismatches
!        Phiall=elliptic3D_pol_conv(srcterm,sig0all,sigPall,sigHall,incapall,v2all,v3all, &
!                          Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
!                          dt,dx1,dx1i,dx2,dx2i,dx3all,dx3iall,Phiall,flagdirich)    !electrodynamic solution
      end if
    else
      write(*,*) 'Skipping potential solution for this time step...'
    end if


    !STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL
    E1prevall=E1all
    E2prevall=E2all
    E3prevall=E3all


    !CALCULATE FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))    !causes major memory leak. maybe from arithmetic statement argument? Left here as a 'lesson learned'
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
    Phiall=-1d0*Phiall
    E1all=grad3D1(Phiall,dx1(1:lx1))
    E2all=grad3D2(Phiall,dx2(1:lx2))
    E3all=grad3D3(Phiall,dx3all(1:lx3all))
    Phiall=-1d0*Phiall   !put things back for good form


    !ADD IN BACKGROUND FIELDS FOR WORKERS
    write(*,*) 'Max field BG and response values are:  ',maxval(pack(abs(E02all),.true.)), &
               maxval(pack(abs(E03all),.true.)),maxval(pack(abs(E2all),.true.)),maxval(pack(abs(E3all),.true.))
    E2all=E2all+E02all
    E3all=E3all+E03all


    !COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
    if (maxval(pack(incap,.true.)) > 1d-6) then
      grad2E=grad3D2(E2all,dx2(1:lx2))
      grad3E=grad3D3(E2all,dx3all(1:lx3all))
      DE2Dtall=(E2all-E2prevall)/dt+v2all*grad2E+v3all*grad3E
      grad2E=grad3D2(E3all,dx2(1:lx2))
      grad3E=grad3D3(E3all,dx3all(1:lx3all))
      DE3Dtall=(E3all-E3prevall)/dt+v2all*grad2E+v3all*grad3E
      J1polall=0d0
      J2polall=incapall*DE2Dtall
      J3polall=incapall*DE3Dtall
    else
      DE2Dtall=0d0
      DE3Dtall=0d0
      J1polall=0d0
      J2polall=0d0
      J3polall=0d0
    end if


    !BROADCAST DERIVED ELECTRIC FIELDS TO THE WORKERS
    call bcast_send(E2all,tagE2,E2)
    call bcast_send(E3all,tagE3,E3)

    !CALCULATE AND DISTRIBUTE CURRENT TO WORKERS
    J1all=sig0all*E1all
    J2all=sigPall*E2all-sigHall*E3all    !BG field already added
    J3all=sigHall*E2all+sigPall*E3all

    write(*,*) 'Max FAC computed to be:  ',maxval(pack(abs(J1all(lx1,:,:)),.true.))
    write(*,*) 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(pack(abs(J2polall),.true.)), &
                 maxval(pack(abs(J3polall),.true.))
    write(*,*) 'Max conduction J2,3 (abs. val.) computed to be:  ',maxval(pack(abs(J2all),.true.)), &
                 maxval(pack(abs(J3all),.true.))


    !GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
    J1all=J1all+J1polall
    J2all=J2all+J2polall
    J3all=J3all+J3polall


    !BROADCAST PARALLEL FIELD AND CURRENT DENSITY TO WORKERS.
    call bcast_send(J1all,tagJ1,J1)
    call bcast_send(E1all,tagE1,E1)

  end subroutine potential_root_mpi_cart


  subroutine potential_root_mpi_curv(it,t,dt,sig0,sigP,sigH,incap,vs2,vs3,x, &
                                potsolve, &
                                E1,E2,E3,J1,J2,J3, &
                                E1all,E2all,E3all,J1all,J2all,J3all,Phiall)

    !------------------------------------------------------------
    !-------ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL
    !------------------------------------------------------------

    integer, intent(in) :: it    !note that the makefile and ilupack insist that the default integer type is 8...
    real(wp), intent(in) :: t,dt
    real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH
    real(wp), dimension(:,:,:), intent(in) ::  incap
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3

    type(curvmesh), intent(in) :: x

    integer, intent(in) :: potsolve

    real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3
    real(wp), dimension(:,:,:), intent(inout) :: E1all,E2all,E3all
    real(wp), dimension(:,:,:), intent(out) :: J1all,J2all,J3all
    real(wp), dimension(:,:,:), intent(inout) :: Phiall   !not good form, but I'm lazy

    real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: sig0all,sigPall,sigHall   !really prefer sizes set by intent(in) vars.
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: incapall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: J1polall,J2polall,J3polall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: E1prevall,E2prevall,E3prevall
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: DE2Dt,DE3Dt

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: srcterm
    real(wp), dimension(1:size(E1all,2),1:size(E1all,3)) :: Vminx1,Vmaxx1
    real(wp), dimension(1:size(E1all,1),1:size(E1all,3)) :: Vminx2,Vmaxx2
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2)) :: Vminx3,Vmaxx3
    integer(8) :: flagdirich    !ilupack expects 8-byte integers with the compile flags presently used

    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: v2all,v3all   !stores drift velocs. for pol. current
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: E01all,E02all,E03all
    real(wp), dimension(1:size(E1all,1),1:size(E1all,2),1:size(E1all,3)) :: DE2Dtall,DE3Dtall   !pol. drift

    integer :: ix2,ix3,lx1,lx2,lx3,lx3all

    real(wp) :: dVmaxx1,dE02all,dE03all

    integer :: iid,islstart,islfin


    !SIZES
    lx1=size(sig0,1)
    lx2=size(sig0,2)
    lx3=size(sig0,3)
    lx3all=size(E1all,3)


    !GATHER WORKER CONDUCTIVITY VALUES
    call gather_recv(sigP,tagsigP,sigPall)
    call gather_recv(sigH,tagsigH,sigHall)
    call gather_recv(sig0,tagsig0,sig0all)
    call gather_recv(incap,tagincap,incapall)
    v2=vs2(1:lx1,1:lx2,1:lx3,1)
    call gather_recv(v2,tagv2pol,v2all)
    v3=vs3(1:lx1,1:lx2,1:lx3,1)
    call gather_recv(v3,tagv3pol,v3all)


    !POPULATE BACKGROUND AND BOUNDARY CONDITION ARRAYS
    call potentialBCs3D(t,sig0all,x%x1,x%x2,x%x3all,x%dx2,x%dx3all,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                        E01all,E02all,E03all,flagdirich)
    J1all=0d0    !so this div is only perp components
    J2all=sigPall*E02all-sigHall*E03all    !BG x2 current
    J3all=sigHall*E02all+sigPall*E03all    !BG x3 current
    srcterm=div3D(J1all,J2all,J3all,x,1,lx1,1,lx2,1,lx3all)


    !FORM AN INITIAL GUESS FOR THE ITERATIVE SOLVER AND CALL SOLVER
    if (it==1) then
      write(*,*) 'Using initial potential guess based on top BC...'
      if (flagdirich==1) then
        Phiall=spread(Vmaxx1,1,lx1)   !very good initial guess
      else
        Phiall=0d0     !will require many iterations to fix
      end if
    else
      write(*,*) 'Using initial potential guess based on previous time step solution...'
    end if


    !SOLVE FOR POTENTIAL USING INCOMPLETE LU FACTORIZATION
    if (mod(it,1)==0 .or. it==1) then    !Don't really need to evaluate every time, usually every 5 steps is okay...
      dVmaxx1=abs(maxval(pack(Vmaxx1,.true.))-minval(pack(Vmaxx1,.true.)))    !gfortran balks at massive if argument...
      dE02all=maxval(pack(abs(E02all),.true.))
      dE03all=maxval(pack(abs(E03all),.true.))

      if (dVmaxx1 < 1d-6 .and. dE02all < 1d-12 .and. dE03all < 1d-12) then   !also no eval if nothing is happening
        write(*,*) 'Potential boundary conditions AND source terms appear null, &
                    moving on with zero potential everywhere...'
        Phiall=0d0
      else
        write(*,*) 'Executing full solve for potential this time step...'
        Phiall=elliptic3D(srcterm,sig0all,sigPall,sigHall,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                        x%dx1,x%dx1i,x%dx2,x%dx2i,x%dx3all,x%dx3iall,Phiall,flagdirich)    !pure electrostatics, note all args are doubles,logicals so no integer type mismatches
!        Phiall=elliptic3D_pol_conv(srcterm,sig0all,sigPall,sigHall,incapall,v2all,v3all, &
!                          Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
!                          dt,dx1,dx1i,dx2,dx2i,dx3all,dx3iall,Phiall,flagdirich)    !electrodynamic solution
      end if
    else
      write(*,*) 'Skipping potential solution for this time step...'
    end if


    !STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL
    E1prevall=E1all
    E2prevall=E2all
    E3prevall=E3all


    !CALCULATE FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))    !causes major memory leak. maybe from arithmetic statement argument? Left here as a 'lesson learned'
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
    Phiall=-1d0*Phiall
    E1all=grad3D1(Phiall,x,1,lx1)
    E2all=grad3D2(Phiall,x,1,lx2)
    E3all=grad3D3(Phiall,x,1,lx3all)
    Phiall=-1d0*Phiall   !put things back for good form


    !ADD IN BACKGROUND FIELDS FOR WORKERS
    write(*,*) 'Max field BG and response values are:  ',maxval(pack(abs(E02all),.true.)), &
               maxval(pack(abs(E03all),.true.)),maxval(pack(abs(E2all),.true.)),maxval(pack(abs(E3all),.true.))
    E2all=E2all+E02all
    E3all=E3all+E03all


    !COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
    if (maxval(pack(incap,.true.)) > 1d-6) then
      grad2E=grad3D2(E2all,x,1,lx2)
      grad3E=grad3D3(E2all,x,1,lx3all)
      DE2Dtall=(E2all-E2prevall)/dt+v2all*grad2E+v3all*grad3E
      grad2E=grad3D2(E3all,x,1,lx2)
      grad3E=grad3D3(E3all,x,1,lx3all)
      DE3Dtall=(E3all-E3prevall)/dt+v2all*grad2E+v3all*grad3E
      J1polall=0d0
      J2polall=incapall*DE2Dtall
      J3polall=incapall*DE3Dtall
    else
      DE2Dtall=0d0
      DE3Dtall=0d0
      J1polall=0d0
      J2polall=0d0
      J3polall=0d0
    end if


    !BROADCAST DERIVED ELECTRIC FIELDS TO THE WORKERS
    call bcast_send(E2all,tagE2,E2)
    call bcast_send(E3all,tagE3,E3)

    !CALCULATE AND DISTRIBUTE CURRENT TO WORKERS
    J1all=sig0all*E1all
    J2all=sigPall*E2all-sigHall*E3all    !BG field already added
    J3all=sigHall*E2all+sigPall*E3all

    write(*,*) 'Max FAC computed to be:  ',maxval(pack(abs(J1all(lx1,:,:)),.true.))
    write(*,*) 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(pack(abs(J2polall),.true.)), &
                 maxval(pack(abs(J3polall),.true.))
    write(*,*) 'Max conduction J2,3 (abs. val.) computed to be:  ',maxval(pack(abs(J2all),.true.)), &
                 maxval(pack(abs(J3all),.true.))


    !GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
    J1all=J1all+J1polall
    J2all=J2all+J2polall
    J3all=J3all+J3polall


    !BROADCAST PARALLEL FIELD AND CURRENT DENSITY TO WORKERS.
    call bcast_send(J1all,tagJ1,J1)
    call bcast_send(E1all,tagE1,E1)

  end subroutine potential_root_mpi_curv

  subroutine potential_workers_mpi(sig0,sigP,sigH,incap,vs2,vs3,potsolve, &
                          E1,E2,E3,J1,J2,J3)

    !------------------------------------------------------------
    !-------WORKER MPI COMMUNICATION ROUTINE FOR POTENTIAL SOLN.
    !------------------------------------------------------------

    real(wp), dimension(:,:,:), intent(in) ::  sig0,sigP,sigH
    real(wp), dimension(:,:,:), intent(in) ::  incap
    real(wp), dimension(-1:,-1:,-1:,:), intent(in) ::  vs2,vs3

    integer, intent(in) :: potsolve

    real(wp), dimension(:,:,:), intent(out) :: E1,E2,E3,J1,J2,J3

    real(wp), dimension(1:size(sigP,1),1:size(sigP,2),1:size(sigP,3)) :: vEB    !send vars need to be contiguous in memory

    integer :: lx1,lx2,lx3

    lx1=size(sig0,1)
    lx2=size(sig0,2)
    lx3=size(sig0,3)


    !SEND DATA TO ROOT, WHO IS GATHERING
    call gather_send(sigP,tagsigP)
    call gather_send(sigH,tagsigH)
    call gather_send(sig0,tagsig0)
    call gather_send(incap,tagincap)
    vEB=vs2(1:lx1,1:lx2,1:lx3,1)    !O+ drift
    call gather_send(vEB,tagv2pol)
    vEB=vs3(1:lx1,1:lx2,1:lx3,1)
    call gather_send(vEB,tagv3pol)


    !RECEIVE RESULTS FROM BROADCASTING ROOT
    call bcast_recv(E2,tagE2)
    call bcast_recv(E3,tagE3)

    call bcast_recv(J1,tagJ1)
    call bcast_recv(E1,tagE1)

  end subroutine potential_workers_mpi


end module potential_comm
