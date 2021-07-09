submodule (potential_comm) potential_worker

implicit none (type, external)

contains

module procedure potential_workers_mpi
  !! ROOT MPI COMM./SOLVE ROUTINE FOR POTENTIAL.  THIS VERSION
  !! INCLUDES THE POLARIZATION CURRENT TIME DERIVATIVE PART
  !! AND CONVECTIVE PARTS IN MATRIX SOLUTION.
  !! STATE VARIABLES VS2,3 INCLUDE GHOST CELLS.  FOR NOW THE
  !! POLARIZATION TERMS ARE PASSED BACK TO MAIN FN, EVEN THOUGH
  !! THEY ARE NOT USED (THEY MAY BE IN THE FUTURE)
  
  !integer :: flagdirich
  
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: paramtrim    !to hold trimmed magnetic field
  
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J1pol,J2pol,J3pol
  
  !real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E01,E02,E03!,E02src,E03src   !distributed background fields
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: srcterm!,divJperp
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E1prev,E2prev,E3prev
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: Phi
  
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: integrand,sigintegral    !general work array for doing integrals
  real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPint2,SigPint3,SigHint,incapint,srctermint
  
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: sig0scaled,sigPscaled,sigHscaled
  
  logical :: perflag    !MUMPS stuff
  
  !real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: Vminx1slab,Vmaxx1slab
  
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3
  real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: v2slab,v3slab
  
  integer :: ix1,ix2,ix3,lx1,lx2,lx3,lx3all, ierr
  integer :: idleft,idright,iddown,idup
  
  real(wp) :: tstart,tfin
  
  integer :: flagsolve

  !! PkI
  real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J2prev
  real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPintstar
  
  
  !SIZES - PERHAPS SHOULD BE TAKEN FROM GRID MODULE INSTEAD OF RECOMPUTED?
  lx1=size(sig0,1)
  lx2=size(sig0,2)
  lx3=size(sig0,3)


  !! PkI:  last time step current should be saved
  J2prev=J2
  
  
  ! this should always be on by default unless the user wants to turn off and recompile; ~10% savings in mumps time *per time step*
  perflag=.true.
  
  
  !call BGfields_boundaries_worker(flagdirich,E01,E02,E03,Vminx1slab,Vmaxx1slab)
  
  
  !> Compute source terms, check Lagrangian flag
  !if (cfg%flaglagrangian) then     ! Lagrangian grid, omit background fields from source terms
  !  E02src=0._wp; E03src=0._wp
  !else                             ! Eulerian grid, use background fields
  !  E02src=E02; E03src=E03
  !end if
  call potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02src,E03src,vn2,vn3,B1,muP,muH,ns,Ts,x, &
                             cfg%flaggravdrift,cfg%flagdiamagnetic,srcterm)
  
  
  !    !ZZZ - DEBUG BY GETTING THE ENTIRE SOURCETERM ARRAY
  !    call gather_send(srcterm,tag%src)
  
  
  !!!!!!!!
  !-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
  !-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
  if (lx2/=1 .and. lx3/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
    if (cfg%potsolve == 1) then    !2D, field-integrated solve
      !-------
      !INTEGRATE CONDUCTANCES AND CAPACITANCES FOR SOLVER COEFFICIENTS
      integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)
      sigintegral=integral3D1(integrand,x,1,lx1)    !no haloing required for a field-line integration
      SigPint2=sigintegral(lx1,:,:)
  
      integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)
      sigintegral=integral3D1(integrand,x,1,lx1)
      SigPint3=sigintegral(lx1,:,:)

      !! PkI
      integrand=sigP*x%h1(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)
      sigintegral=integral3D1(integrand,x,1,lx1)
      SigPintstar=sigintegral(lx1,:,:)
  
      integrand=x%h1(1:lx1,1:lx2,1:lx3)*sigH
      sigintegral=integral3D1(integrand,x,1,lx1)
      SigHint=sigintegral(lx1,:,:)
  
      sigintegral=integral3D1(incap,x,1,lx1)
      incapint=sigintegral(lx1,:,:)
      !-------
  
      !PRODUCE A FIELD-INTEGRATED SOURCE TERM
      if (flagdirich /= 1) then
        !! Neumann conditions; incorporate a source term and execute the solve
        !-------
        integrand = x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)*srcterm
        sigintegral = integral3D1(integrand,x,1,lx1)
        srctermint = sigintegral(lx1,:,:)
        srctermint = srctermint+x%h2(lx1,1:lx2,1:lx3)*x%h3(lx1,1:lx2,1:lx3)*Vmaxx1slab- &
                                   x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab
        !! workers don't have access to boundary conditions, unless root sends
        !-------
  
  
        !RADD--- ROOT NEEDS TO PICK UP *INTEGRATED* SOURCE TERMS AND COEFFICIENTS FROM WORKERS
        call gather_send(srctermint,tag%src)
        call gather_send(incapint,tag%incapint)
        call gather_send(SigPint2,tag%SigPint2)
        call gather_send(SigPint3,tag%SigPint3)
        call gather_send(SigHint,tag%SigHint)
        v2=vs2(1:lx1,1:lx2,1:lx3,1); v3=vs3(1:lx1,1:lx2,1:lx3,1);
        v2slab=v2(lx1,:,:); v3slab=v3(lx1,:,:)
        call gather_send(v2slab,tag%v2electro)
        call gather_send(v3slab,tag%v3electro)
  !          v2slab=vs2(lx1,1:lx2,1:lx3,1); v3slab=vs3(lx1,1:lx2,1:lx3,1);
        !! need to pick out the ExB drift here (i.e. the drifts from highest altitudes);
        !! but this is only valid for Cartesian, so it's okay for the foreseeable future

        !! PkI
        call gather_send(SigPintstar,tagSigPint2)
        call gather_send(J2prev,tagJ2)
  
  
        call elliptic_workers()    !workers do not need any specific info about the problem (that all resides with root who will redistribute)
      else
        !! Dirichlet conditions
        !! - since this is field integrated we just copy BCs specified by user to other locations along field line (root does this)
  
      end if
  
  !
    else    !resolved 3D solve
      !! ZZZ - conductivities need to be properly scaled here...
      !! So does the source term...  Maybe leave as broken for now since I don't really plan to use this code
     !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
      sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
      srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
      sigHscaled=x%h1(1:lx1,1:lx2,1:lx3)*sigH
  
      !RADD--- ROOT NEEDS TO PICK UP FIELD-RESOLVED SOURCE TERM AND COEFFICIENTS FROM WORKERS
      call gather_send(sigPscaled,tag%sigP)
      call gather_send(sigHscaled,tag%sigH)
      call gather_send(sig0scaled,tag%sig0)
      call gather_send(srcterm,tag%src)
  
      call mpi_recv(flagsolve,1,MPI_INTEGER,0,tag%flagdirich,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  
      if (flagsolve/=0) then
        call elliptic_workers()
      end if
  
    end if
  else   !lx1=1 so do a field-resolved 2D solve over x1,x3
  
  
    !-------
    !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
    sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
    srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
    !-------
  
    !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
    call gather_send(sigPscaled,tag%sigP)
    call gather_send(sig0scaled,tag%sig0)
    call gather_send(srcterm,tag%src)
  
    ! Need to convert current boundary condition into potential normal derivative
    if (flagdirich==0) then
      if (gridflag==1) then
        Vminx1slab=-x%h1(1,1:lx2,1:lx3)*Vminx1slab/sig0(1,:,:)
        call gather_send(Vminx1slab,tag%Vminx1)
      else
        Vmaxx1slab=-x%h1(lx1,1:lx2,1:lx3)*Vmaxx1slab/sig0(lx1,:,:)
        call gather_send(Vmaxx1slab,tag%Vmaxx1)
      end if
    end if
  
    call elliptic_workers()
  end if
  !    print *, 'MUMPS time:  ',tfin-tstart
  !!!!!!!!!
  
  
  !RADD--- ROOT NEEDS TO PUSH THE POTENTIAL BACK TO ALL WORKERS FOR FURTHER PROCESSING (BELOW)
  call bcast_recv(Phi,tag%Phi)
  
  
  !-------
  !! STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL
  !! (OLD FIELDS USED TO CALCULATE POLARIZATION CURRENT)
  E1prev=E1
  E2prev=E2
  E3prev=E3
  !-------
  
  
  !-------
  !CALCULATE PERP FIELDS FROM POTENTIAL
  !      E20all=grad3D2(-Phi0all,dx2(1:lx2))
  !! causes major memory leak. maybe from arithmetic statement argument?
  !! Left here as a 'lesson learned' (or is it a gfortran bug...)
  !      E30all=grad3D3(-Phi0all,dx3all(1:lx3all))
  call pot2perpfield(Phi,x,E2,E3)
  !--------
  
  
  !    !R-------
  !    !JUST TO JUDGE THE IMPACT OF MI COUPLING
  !    print *, 'Max integrated inertial capacitance:  ',maxval(incapint)
  !    print *, 'Max integrated Pedersen conductance (includes metric factors):  ',maxval(SigPint2)
  !    print *, 'Max integrated Hall conductance (includes metric factors):  ',minval(SigHint), maxval(SigHint)
  !!    print *, 'Max E2,3 BG and response values are:  ',maxval(abs(E02)), maxval(abs(E03)), maxval(abs(E2)),maxval(abs(E3))
  !    print *, 'Max E2,3 BG and response values are:  ',maxval(E02), maxval(E03),maxval(E2),maxval(E3)
  !    print *, 'Min E2,3 BG and response values are:  ',minval(E02), minval(E03),minval(E2),minval(E3)
  !    !R-------
  
  
  !--------
  !ADD IN BACKGROUND FIELDS BEFORE HALOING
  !if (.not. cfg%flaglagrangian) then
  !  E2=E2+E02
  !  E3=E3+E03
  !end if
  !--------
  
  
  call polarization_currents(cfg,x,dt,incap,E2,E3,E2prev,E3prev,v2,v3,J1pol,J2pol,J3pol)
  
  
  !--------
  J2=0._wp; J3=0._wp    ! must be zeroed out before we accumulate currents
  call acc_perpconductioncurrents(sigP,sigH,E2,E3,J2,J3)
  call acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)
  if (cfg%flagdiamagnetic) then
    call acc_pressurecurrents(muP,muH,ns,Ts,x,J2,J3)
  end if
  if (cfg%flaggravdrift) then
    call acc_perpgravcurrents(sigPgrav,sigHgrav,g2,g3,J2,J3)
  end if
  !--------
  
  
  call parallel_currents(cfg,x,J2,J3,Vminx1slab,Vmaxx1slab,Phi,sig0,flagdirich,J1,E1)
  
  
  !    !R-------
  if (debug) then
  !    print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,:,:)))
  !! ZZZ - this rey needsz to be current at the "top"
  !    print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)),  maxval(abs(J3pol))
  !    print *, 'Max conduction J2,3  computed to be:  ',maxval(J2), maxval(J3)
  !    print *, 'Min conduction J2,3  computed to be:  ',minval(J2),  minval(J3)
  !    print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1))
  endif
  !    !R-------
  
  
  !-------
  !GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
  J1=J1+J1pol
  J2=J2+J2pol
  J3=J3+J3pol
  !-------
end procedure potential_workers_mpi

end submodule potential_worker
