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
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: J1pol,J2pol,J3pol
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: srcterm!,divJperp
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: E1prev,E2prev,E3prev
  real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: Phi    ! FIXME: why a local copy???
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: integrand,sigintegral    !general work array for doing integrals
  real(wp), dimension(1:lx2,1:lx3) :: SigPint2,SigPint3,SigHint,incapint,srctermint
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: sig0scaled,sigPscaled,sigHscaled
  logical :: perflag    !MUMPS stuff
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: v2,v3
  real(wp), dimension(1:lx2,1:lx3) :: v2slab,v3slab
  integer :: flagsolve

  ! this should always be on by default unless the user wants to turn off and recompile; ~10% savings in mumps time *per time step*
  perflag=.false.
  call potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02src,E03src,vn2,vn3,B1,muP,muH,ns,Ts,x, &
                             cfg%flaggravdrift,cfg%flagdiamagnetic,cfg%flagnodivJ0,srcterm)

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

      call mpi_recv(flagsolve,1,MPI_INTEGER,0,tag%flagdirich,MPI_COMM_WORLD,MPI_STATUS_IGNORE)

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
  call bcast_recv3D_ghost(Phi,tag%Phi)

  !-------
  !! STORE PREVIOUS TIME TOTAL FIELDS BEFORE UPDATING THE ELECTRIC FIELDS WITH NEW POTENTIAL
  !! (OLD FIELDS USED TO CALCULATE POLARIZATION CURRENT)
  E1prev=E1(1:lx1,1:lx2,1:lx3)
  E2prev=E2(1:lx1,1:lx2,1:lx3)
  E3prev=E3(1:lx1,1:lx2,1:lx3)
  !-------

  !-------
  !CALCULATE PERP FIELDS FROM POTENTIAL
  !      E20all=grad3D2(-Phi0all,dx2(1:lx2))
  !! causes major memory leak. maybe from arithmetic statement argument?
  !! Left here as a 'lesson learned' (or is it a gfortran bug...)
  !      E30all=grad3D3(-Phi0all,dx3all(1:lx3all))
  call pot2perpfield(Phi,x,E2,E3)
  !--------

  call polarization_currents(cfg,x,dt,incap,E2,E3,E2prev,E3prev,v2,v3,J1pol,J2pol,J3pol)

  !--------
  J2(1:lx1,1:lx2,1:lx3)=0._wp; J3(1:lx1,1:lx2,1:lx3)=0._wp    ! must be zeroed out before we accumulate currents
  if (.not. cfg%flagnodivJ0) call acc_perpBGconductioncurrents(sigP,sigH,E02src,E03src,J2,J3)
  !^ note that out input background fields to this procedure have already been tweaked to account for lagrangian vs. eulerian grids so we can just blindly add these in without worry
  call acc_perpconductioncurrents(sigP,sigH,E2,E3,J2,J3)
  call acc_perpwindcurrents(sigP,sigH,vn2,vn3,B1,J2,J3)
  if (cfg%flagdiamagnetic) then
    call acc_pressurecurrents(muP,muH,ns,Ts,x,J2,J3)
  end if
  if (cfg%flaggravdrift) then
    call acc_perpgravcurrents(sigPgrav,sigHgrav,x%g2,x%g3,J2,J3)
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
  J1(1:lx1,1:lx2,1:lx3)=J1(1:lx1,1:lx2,1:lx3)+J1pol
  J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+J2pol
  J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)+J3pol
  !-------
end procedure potential_workers_mpi

end submodule potential_worker
