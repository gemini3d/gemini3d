submodule (potential_comm) potential_root

implicit none (type, external)


contains


module procedure potential_root_mpi_curv

!! ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL.  THIS VERSION
!! INCLUDES THE POLARIZATION CURRENT TIME DERIVATIVE PART
!! AND CONVECTIVE PARTS IN MATRIX SOLUTION.
!! STATE VARIABLES VS2,3 INCLUDE GHOST CELLS.  FOR NOW THE
!! POLARIZATION TERMS ARE PASSED BACK TO MAIN FN, EVEN THOUGH
!! THEY ARE NOT USED (THEY MAY BE IN THE FUTURE)

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3

real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: srctermall
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)), target :: Vminx1,Vmaxx1     !allow pointer aliases for these vars.
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: Vminx1buf,Vmaxx1buf
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,3)) :: Vminx2,Vmaxx2
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2)) :: Vminx3,Vmaxx3
integer :: flagdirich

real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: v2slaball,v3slaball   !stores drift velocs. for pol. current

!   real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: paramtrim    !to hold trimmed magnetic field
real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: E01all,E02all,E03all    !background fields
!   real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: divJperpall2,divJperpall3,divJperpall
!! more work arrays

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: integrand,sigintegral    !general work array for doing integrals

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: DE2Dt,DE3Dt   !pol. drift
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J1pol,J2pol,J3pol

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E01,E02,E03   !distributed background fields
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: srcterm,divJperp
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E1prev,E2prev,E3prev
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: Phi

real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPint2,SigPint3,SigHint,incapint,srctermint
real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: SigPint2all,SigPint3all,SigHintall,incapintall,srctermintall

real(wp), dimension(1:size(Phiall,2),1:size(Phiall,3)) :: Phislab,Phislab0

real(wp), dimension(0:size(E1,1)+1,0:size(E1,2)+1,0:size(E1,3)+1) :: divtmp
!! one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(E1,1)+2,-1:size(E1,2)+2,-1:size(E1,3)+2) :: J1halo,J2halo,J3halo
!! haloing assumes existence of two ghost cells

real(wp), dimension(1:size(Phiall,1),1:size(Phiall,2),1:size(Phiall,3)) :: sig0scaledall,sigPscaledall,sigHscaledall
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: sig0scaled,sigPscaled,sigHscaled

logical :: perflag    !MUMPS stuff

real(wp), dimension(1:size(Phiall,3)) :: Vminx2slice,Vmaxx2slice
real(wp), dimension(1:size(Phiall,2)) :: Vminx3slice,Vmaxx3slice
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: Vminx1slab,Vmaxx1slab
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: v2slab,v3slab

integer :: iid, ierr
integer :: ix1,ix2,ix3,lx1,lx2,lx3,lx3all
integer :: idleft,idright,iddown,idup

real(wp) :: tstart,tfin


!SIZES - PERHAPS SHOULD BE TAKEN FROM GRID MODULE INSTEAD OF RECOMPUTED?
lx1=size(sig0,1)
lx2=size(sig0,2)
lx3=size(sig0,3)
lx3all=size(Phiall,3)


!USE PREVIOUS MUMPS PERMUTATION (OLD CODE? BUT MIGHT BE WORTH REINSTATING?)
perflag=.true.


!R-------
!! POPULATE BACKGROUND AND BOUNDARY CONDITION ARRAYS
!! - IDEALLY ROOT ONLY SINCE IT INVOLVES FILE INPUT, since the interpolation is 2D don't distribute to workers...
call cpu_time(tstart)
if (cfg%flagE0file==1) then
  call potentialBCs2D_fileinput(dt,t,ymd,UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      E01all,E02all,E03all,flagdirich)
else
  call potentialBCs2D(UTsec,cfg,x,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                      E01all,E02all,E03all,flagdirich)     !user needs to manually swap x2 and x3 in this function, e.g. for EIA, etc.
end if
call cpu_time(tfin)
if (debug) print *, 'Root has computed BCs in time:  ',tfin-tstart
!R-------


!R--------
ierr=0
do iid=1,lid-1    !communicate intent for solve to workers so they know whether or not to call mumps fn.
  call mpi_send(flagdirich,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD,ierr)
end do
if (ierr/=0) error stop 'mpi_send failed to send solve intent'
if (debug) print *, 'Root has communicated type of solve to workers:  ',flagdirich
!R--------


!Need to broadcast background fields from root
!Need to also broadcast x1 boundary conditions for source term calculations.
call bcast_send(E01all,tag%E01,E01)
call bcast_send(E02all,tag%E02,E02)
call bcast_send(E03all,tag%E03,E03)

!These are pointer targets so don't assume contiguous in memory - pack them into a buffer to be safe
Vmaxx1buf=Vmaxx1; Vminx1buf=Vminx1;
call bcast_send(Vminx1buf,tag%Vminx1,Vminx1slab)
call bcast_send(Vmaxx1buf,tag%Vmaxx1,Vmaxx1slab)


!> Compute source terms
call potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02,E03,vn2,vn3,B1,x,cfg%flaggravdrift,srcterm)


!!!!!!!!
!-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
!-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
if (lx2/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
  if (cfg%potsolve == 1) then    !2D, field-integrated solve
    if (debug) print *, 'Beginning field-integrated solve...'


    !> INTEGRATE CONDUCTANCES AND CAPACITANCES FOR SOLVER COEFFICIENTS
    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)    !no haloing required for a field-line integration
    SigPint2=sigintegral(lx1,:,:)

    integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigPint3=sigintegral(lx1,:,:)

    integrand=x%h1(1:lx1,1:lx2,1:lx3)*sigH
    sigintegral=integral3D1(integrand,x,1,lx1)
    SigHint=sigintegral(lx1,:,:)

    sigintegral=integral3D1(incap,x,1,lx1)
    incapint=sigintegral(lx1,:,:)
    !-------


    !> PRODUCE A FIELD-INTEGRATED SOURCE TERM
    if (flagdirich /= 1) then   !Neumann conditions; incorporate a source term and execute the solve
      if (debug) print *, 'Using FAC boundary condition...'
      !-------
      integrand=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*srcterm
      sigintegral=integral3D1(integrand,x,1,lx1)
      srctermint=sigintegral(lx1,:,:)
      srctermint=srctermint+x%h2(lx1,1:lx2,1:lx3)*x%h3(lx1,1:lx2,1:lx3)*Vmaxx1slab- &
                                 x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab
      !! workers don't have access to boundary conditions, unless root sends
      !-------


      v2=vs2(1:lx1,1:lx2,1:lx3,1); v3=vs3(1:lx1,1:lx2,1:lx3,1);
      ! must be set since used later by the polarization current calculation
      v2slab=v2(lx1,1:lx2,1:lx3); v3slab=v3(lx1,1:lx2,1:lx3);
      !! need to pick out the ExB drift here (i.e. the drifts from highest altitudes); but this is only valid for Cartesian,
      !! so it's okay for the foreseeable future


      !RADD--- ROOT NEEDS TO PICK UP *INTEGRATED* SOURCE TERMS AND COEFFICIENTS FROM WORKERS
      call gather_recv(srctermint,tag%src,srctermintall)
      call gather_recv(incapint,tag%incapint,incapintall)
      call gather_recv(SigPint2,tag%SigPint2,SigPint2all)
      call gather_recv(SigPint3,tag%SigPint3,SigPint3all)
      call gather_recv(SigHint,tag%SigHint,SigHintall)
      call gather_recv(v2slab,tag%v2electro,v2slaball)
      call gather_recv(v3slab,tag%v3electro,v3slaball)


      !R------
      !EXECUTE FIELD-INTEGRATED SOLVE
      Vminx2slice=Vminx2(lx1,:)    !slice the boundaries into expected shape
      Vmaxx2slice=Vmaxx2(lx1,:)
      Vminx3slice=Vminx3(lx1,:)
      Vmaxx3slice=Vmaxx3(lx1,:)
      Phislab0=Phiall(lx1,:,:)    !root already possess the fullgrid potential from prior solves...
      if (debug) print *, 'Root is calling MUMPS...'
      !R-------

      !R------ EXECUTE THE MUMPS SOLVE FOR FIELD-INT
      call cpu_time(tstart)
      if (.not. x%flagper) then     !nonperiodic mesh
        if (debug) print *, '!!!User selected aperiodic solve...'
        Phislab=potential2D_polarization(srctermintall,SigPint2all,SigPint3all,SigHintall,incapintall,v2slaball,v3slaball, &
                                 Vminx2slice,Vmaxx2slice,Vminx3slice,Vmaxx3slice, &
                                 dt,x,Phislab0,perflag,it)
        !! note tha this solver is only valid for cartesian meshes, unless the inertial capacitance is set to zero
      else
        if (debug) print *, '!!!User selected periodic solve...'
        Phislab = potential2D_polarization_periodic(srctermintall,SigPint2all,SigHintall,incapintall,v2slaball,v3slaball, &
                                   Vminx2slice,Vmaxx2slice,Vminx3slice,Vmaxx3slice, &
                                   dt,x,Phislab0,perflag,it)
        !! !note that either sigPint2 or 3 will work since this must be cartesian...
      end if
      call cpu_time(tfin)
      if (debug) print *, 'Root received results from MUMPS which took time:  ',tfin-tstart
      !R-------

    else
      !! Dirichlet conditions - since this is field integrated we just copy BCs specified by user
      !! to other locations along field line (done later)
      !R------
      Phislab=Vmaxx1
      !! potential is whatever user specifies, since we assume equipotential field lines,
      !! it doesn't really matter whether we use Vmaxx1 or Vminx1.
      !! Note however, tha thte boundary conditions subroutines should explicitly
      !! set these to be equal with Dirichlet conditions, for consistency.
      if (debug) print *, 'Dirichlet conditions selected with field-integrated solve. Copying BCs along x1-direction...'
      !R------
    end if


    !R------  AFTER ANY TYPE OF FIELD-INT SOLVE COPY THE BCS ACROSS X1 DIMENSION
    do ix1=1,lx1
      Phiall(ix1,:,:)=Phislab(:,:)
      !! copy the potential across the ix1 direction; past this point there is no real difference with 3D,
      !! note that this is still valid in curvilinear form
    end do
    !R------

  else
    !! resolved 3D solve
    !! ZZZ - conductivities need to be properly scaled here...
    !! So does the source term...  Maybe leave as broken for now since there are no immediate plans to use this (too slow)
    if (debug) print *, 'Beginning field-resolved 3D solve...  Type;  ',flagdirich

    !-------
    !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM,
    !need to adopt for curvilinear case...
    sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
    if (flagswap==1) then
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
    else
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
    end if
    srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
    sigHscaled=x%h1(1:lx1,1:lx2,1:lx3)*sigH


    !RADD--- ROOT NEEDS TO PICK UP FIELD-RESOLVED SOURCE TERM AND COEFFICIENTS FROM WORKERS
    call gather_recv(sigPscaled,tag%sigP,sigPscaledall)
    call gather_recv(sigHscaled,tag%sigH,sigHscaledall)
    call gather_recv(sig0scaled,tag%sig0,sig0scaledall)
    call gather_recv(srcterm,tag%src,srctermall)


    !R------
    if (debug) print *, '!Beginning field-resolved 3D solve (could take a very long time)...'
   ! Phiall=elliptic3D_curv(srctermall,sig0scaledall,sigPscaledall,sigHscaledall,Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
   !                   x,flagdirich,perflag,it)
    if( maxval(abs(Vminx1))>1e-12_wp .or. maxval(abs(Vmaxx1))>1e-12_wp ) then
      do iid=1,lid-1
        call mpi_send(1,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD,ierr)
      end do
      Phiall=potential3D_fieldresolved_decimate(srctermall,sig0scaledall,sigPscaledall,sigHscaledall, &
                                 Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                 x,flagdirich,perflag,it)
    else
      do iid=1,lid-1
        call mpi_send(0,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD,ierr)
      end do
        if (debug) print*, 'Boundary conditions too small to require solve, setting everything to zero...'
      Phiall=0e0_wp
    end if
    !R------
  end if
else   !lx1=1 so do a field-resolved 2D solve over x1,x3
  if (debug) print *, 'Beginning field-resolved 2D solve...  Type;  ',flagdirich


  !-------
  !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
  sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
  if (flagswap==1) then
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
  else
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
  end if
  srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
!      srcterm=-1d0*srcterm
  !! in a 2D solve negate this due to it being a cross produce and the fact that we've permuted the 2 and 3 dimensions.
  !! ZZZ - NOT JUST THIS WORKS WITH BACKGROUND FIELDS???
  !-------

  !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
  call gather_recv(sigPscaled,tag%sigP,sigPscaledall)
  call gather_recv(sig0scaled,tag%sig0,sig0scaledall)
  call gather_recv(srcterm,tag%src,srctermall)


  !! EXECUTE THE SOLVE WITH MUMPS AND SCALED TERMS
  !! NOTE THE LACK OF A SPECIAL CASE HERE TO CHANGE THE POTENTIAL PROBLEM
  !! - ONLY THE HALL TERM CHANGES (SINCE RELATED TO EXB) BUT THAT DOESN'T APPEAR IN THIS EQN!
  Phiall=potential2D_fieldresolved(srctermall,sig0scaledall,sigPscaledall,Vminx1,Vmaxx1,Vminx3,Vmaxx3, &
                    x,flagdirich,perflag,it)
end if
if (debug) print *, 'MUMPS time:  ',tfin-tstart
!!!!!!!!!


!RADD--- ROOT NEEDS TO PUSH THE POTENTIAL BACK TO ALL WORKERS FOR FURTHER PROCESSING (BELOW)
!Phiall=0._wp
call bcast_send(Phiall,tag%Phi,Phi)


!-------
!! STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL
!! (OLD FIELDS USED TO CALCULATE POLARIZATION CURRENT)
E1prev=E1
E2prev=E2
E3prev=E3
!-------


!-------
!CALCULATE PERP FIELDS FROM POTENTIAL
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))
!! causes major memory leak. maybe from arithmetic statement argument?
!! Left here as a 'lesson learned' (or is it a gfortran bug...)
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
call pot2perpfield(Phi,x,E2,E3)


!R-------
!JUST TO JUDGE THE IMPACT OF MI COUPLING
if (debug) then
print *, 'Max integrated inertial capacitance:  ',maxval(incapintall)
!print *, 'Max integrated Pedersen conductance (includes metric factors):  ',maxval(SigPint2all)
!print *, 'Max integrated Hall conductance (includes metric factors):  ',minval(SigHintall), maxval(SigHintall)
print *, 'Max E2,3 BG and response values are:  ',maxval(E02), maxval(E03),maxval(E2),maxval(E3)
print *, 'Min E2,3 BG and response values are:  ',minval(E02), minval(E03),minval(E2),minval(E3)
print *, 'Min/Max values of potential:  ',minval(Phi),maxval(Phi)
print *, 'Min/Max values of full grid potential:  ',minval(Phiall),maxval(Phiall)
endif
!R-------


!--------
!ADD IN BACKGROUND FIELDS BEFORE HALOING
E2=E2+E02
E3=E3+E03
!--------


!--------
!COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
!if (maxval(incap) > 0._wp) then
!! ZZZ this is really bad needs to be a global test rather than having each worker test since there
!! is message passing embedded in here and everyone needs to do the same thing!!!
if (cfg%flagcap/=0) then
  if (debug) print*, 'Working on polarization currents...'

  !differentiate E2 in x2 (needs haloing)
  J1halo(1:lx1,1:lx2,1:lx3)=E2
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D2(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad2E=divtmp(1:lx1,1:lx2,1:lx3)

  !Now differentiate E2 in the x3 direction
  J1halo(1:lx1,1:lx2,1:lx3)=E2
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  !total derivative needed to compute the x2 component of the polarization current
  DE2Dt=(E2-E2prev)/dt+v2*grad2E+v3*grad3E

  !differentiate E3 in x2
  J1halo(1:lx1,1:lx2,1:lx3)=E3
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D2(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad2E=divtmp(1:lx1,1:lx2,1:lx3)

  !differentiate E3 in x3
  J1halo(1:lx1,1:lx2,1:lx3)=E3
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  !total derivative needed for x3 component of pol. current
  DE3Dt=(E3-E3prev)/dt+v2*grad2E+v3*grad3E

  !convert derivative to current density
  J1pol=0d0
  J2pol=incap*DE2Dt
  J3pol=incap*DE3Dt
else       !pure electrostatic solve was done
  if (debug) print*, 'Electrostatics used, skipping polarization current...'
  DE2Dt=0d0
  DE3Dt=0d0
  J1pol=0d0
  J2pol=0d0
  J3pol=0d0
end if
!--------


!--------
J2=0._wp; J3=0._wp    ! must be zeroed out before we accumulate currents
call acc_perpconductioncurrents(sigP,sigH,E2,E3,J2,J3)
call acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)
if (cfg%flaggravdrift) then
  call acc_perpgravcurrents(sigPgrav,sigHgrav,g2,g3,J2,J3)
end if
!--------


!!!!!!!!
!NOW DEAL WITH THE PARALLEL FIELDS AND ALL CURRENTS
if (lx2/=1 .and. cfg%potsolve ==1) then    !we did a field-integrated solve above
  if (debug) print*, 'Appear to need to differentiate to get J1...'

  !-------
  !NOTE THAT A DIRECT E1ALL CALCULATION WILL GIVE ZERO, SO USE INDIRECT METHOD, AS FOLLOWS
  J1=0d0    !a placeholder so that only the perp divergence is calculated - will get overwritten later.
!      divJperp=div3D(J1,J2,J3,x,1,lx1,1,lx2,1,lx3)

  if (cfg%flagJpar) then   ! user can elect not to compute Jpar, which can be prone to artifacts particularly at low resolution
    J1halo(1:lx1,1:lx2,1:lx3)=J1
    J2halo(1:lx1,1:lx2,1:lx3)=J2
    J3halo(1:lx1,1:lx2,1:lx3)=J3

    call halo_pot(J1halo,tag%J1,x%flagper,.false.)
    call halo_pot(J2halo,tag%J2,x%flagper,.false.)
    call halo_pot(J3halo,tag%J3,x%flagper,.false.)

    divtmp=div3D(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),J2halo(0:lx1+1,0:lx2+1,0:lx3+1), &
                 J3halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
    divJperp=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*divtmp(1:lx1,1:lx2,1:lx3)
    if (flagdirich /= 1) then
      !! Neumann conditions, this is boundary location-agnostic since both bottom and top FACs are known
      !! - they have to  be loaded into VVmaxx1 and Vminx1
      if (debug) print *, 'Nuemann boundaries, integrated from highest altitude down to preserve accuracy...'
      if (gridflag==0) then     !closed dipole grid, really would be best off integrating from the source hemisphere
        if (debug) print *,  'Closed dipole grid; integration starting in source hemisphere (if applicable)...', &
                       minval(Vmaxx1slab), &
                       maxval(Vmaxx1slab)
        if (cfg%sourcemlat>=0d0) then    !integrate from northern hemisphere
          if (debug) print *, 'Source is in northern hemisphere (or there is no source)...'
          J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
          do ix1=1,lx1
            J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                             (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
          end do
        else
          if (debug) print *, 'Source in southern hemisphere...'
          J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current starting from minx1
          do ix1=1,lx1
            J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                             (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
          end do
        end if
      elseif (gridflag==1) then    !this would be an inverted grid, this max altitude corresponds to the min value of x1
        if (debug) print *,  'Inverted grid; integration starting at min x1 (highest alt. or southern hemisphere)...', &
                       minval(Vminx1slab), &
                       maxval(Vminx1slab)
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
        end do
      else        !minx1 is at teh bottom of the grid to integrate from max x1
        if (debug) print *,  'Non-inverted grid; integration starting at max x1...', minval(Vmaxx1slab), maxval(Vmaxx1slab)
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
        end do
      end if
    else
      !! Dirichlet conditions - we need to integrate from the ***lowest altitude***
      !! (where FAC is known to be zero, note this is not necessarilty the logical bottom of the grid), upwards (to where it isn't)
      if (gridflag/=2) then    !inverted grid (logical top is the lowest altitude)
        if (debug) print *, 'Inverted grid detected - integrating logical top downward to compute FAC...'
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (J1(ix1,:,:))    !FAC AT TOP ASSUMED TO BE ZERO
        end do
      else      !non-inverted grid (logical bottom is the lowest altitude - so integrate normy)
        if (debug) print *, 'Non-inverted grid detected - integrating logical bottom to top to compute FAC...'
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (-1d0*J1(ix1,:,:))    !FAC AT THE BOTTOM ASSUMED TO BE ZERO
        end do
      end if
    end if
    E1=J1/sig0
    !-------
  end if ! flagJpar
else   !we resolved the field line (either 2D solve or full 3D) so just differentiate normally

  !-------
  Phi=-1d0*Phi
  E1=grad3D1(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required since x1-derivative
  Phi=-1d0*Phi
  J1=sig0*E1
  !-------

end if
!!!!!!!!!


!R-------
if (debug) then
  print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,:,:)))    !ZZZ - this rey needsz to be current at the "top"
  print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)), &
               maxval(abs(J3pol))
  !    print *, 'Max conduction J2,3 (abs. val.) computed to be:  ',maxval(abs(J2)), &
  !                 maxval(abs(J3))
  print *, 'Max conduction J2,3  computed to be:  ',maxval(J2), &
               maxval(J3)
  print *, 'Min conduction J2,3  computed to be:  ',minval(J2), &
               minval(J3)
  print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1))
  print *, 'flagswap:  ',flagswap
endif
!R-------


!-------
!GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
J1=J1+J1pol
J2=J2+J2pol
J3=J3+J3pol
!-------

! if (t>11) then
!   open(newunit=utrace, form='unformatted', access='stream',file='Phiall.raw8', status='replace', action='write')
!   write(utrace) Phiall
!   close(utrace)
!   error stop 'DEBUG'
! end if

end procedure potential_root_mpi_curv

end submodule potential_root
