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
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: v2,v3
  real(wp), dimension(1:lx1,1:lx2all,1:lx3all) :: srctermall
  real(wp), dimension(1:lx2all,1:lx3all) :: v2slaball,v3slaball   !stores drift velocs. for pol. current
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: integrand,sigintegral    !general work array for doing integrals
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: J1pol,J2pol,J3pol
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: srcterm!,divJperp
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: E1prev,E2prev,E3prev
  real(wp), dimension(-1:lx1+2,-1:lx2+2,-1:lx3+2) :: Phi    ! FIXME: why is this a local copy? main program doesn't need?
  real(wp), dimension(1:lx2,1:lx3) :: SigPint2,SigPint3,SigHint,incapint,srctermint
  real(wp), dimension(1:lx2all,1:lx3all) :: SigPint2all,SigPint3all,SigHintall,incapintall,srctermintall
  real(wp), dimension(1:lx2all,1:lx3all) :: Phislab,Phislab0
  real(wp), dimension(1:lx1,1:lx2all,1:lx3all) :: sig0scaledall,sigPscaledall,sigHscaledall
  real(wp), dimension(1:lx1,1:lx2,1:lx3) :: sig0scaled,sigPscaled,sigHscaled
  logical :: perflag    !MUMPS stuff
  real(wp), dimension(1:lx3all) :: Vminx2slice,Vmaxx2slice
  real(wp), dimension(1:lx2all) :: Vminx3slice,Vmaxx3slice
  !real(wp), dimension(1:lx2,1:lx3) :: Vminx1slab,Vmaxx1slab
  real(wp), dimension(1:lx2,1:lx3) :: v2slab,v3slab
  real(wp), dimension(1:lx1,1:lx2all,1:lx3all) :: Phitmp
  integer :: iid
  integer :: ix1
  real(wp) :: tstart,tfin

  !> store a cached ordering for later use (improves performance substantially)
  perflag=.true.

  call potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02src,E03src,vn2,vn3,B1,muP,muH,ns,Ts,x, &
                             cfg%flaggravdrift,cfg%flagdiamagnetic,cfg%flagnodivJ0,srcterm)


  !!!!!!!!
  !-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
  !-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
  if (lx2all/=1 .and. lx3all/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
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
        Phislab0=Phiall(lx1,1:lx2all,1:lx3all)    !root already possess the fullgrid potential from prior solves...
        if (debug) print *, 'Root is calling MUMPS...'
        !R-------

        !R------ EXECUTE THE MUMPS SOLVE FOR FIELD-INT
        call cpu_time(tstart)
        if (.not. x%flagper) then     !nonperiodic mesh
          if (debug) print *, '!!!User selected aperiodic solve...'
          Phislab=potential2D_polarization(srctermintall,SigPint2all,SigPint3all, &
                                   SigHintall,incapintall,v2slaball,v3slaball, &
                                   Vminx2slice,Vmaxx2slice,Vminx3slice,Vmaxx3slice, &
                                   dt,x,Phislab0,perflag,it)
          !! note that this solver is only valid for cartesian meshes, unless the inertial capacitance is set to zero
        else
          if (debug) print *, '!!!User selected periodic solve...'
          Phislab = potential2D_polarization_periodic(srctermintall,SigPint2all,SigHintall, &
                                     incapintall,v2slaball,v3slaball, &
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
        !! Note however, that the boundary conditions subroutines should explicitly
        !! set these to be equal with Dirichlet conditions, for consistency.
        if (debug) print *, 'Dirichlet conditions selected with field-integrated solve. Copying BCs along x1-direction...'
        !R------
      end if

      !R------  AFTER ANY TYPE OF FIELD-INT SOLVE COPY THE BCS ACROSS X1 DIMENSION
      do ix1=1,lx1
        Phiall(ix1,1:lx2all,1:lx3all)=Phislab(1:lx2all,1:lx3all)
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
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
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
      !if( maxval(abs(Vminx1))>1e-12_wp .or. maxval(abs(Vmaxx1))>1e-12_wp ) then
        do iid=1,mpi_cfg%lid-1
          call mpi_send(1,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD)
        end do
        !Phiall=potential3D_fieldresolved_decimate(srctermall,sig0scaledall,sigPscaledall,sigHscaledall, &
        !                           Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
        !                           x,flagdirich,perflag,it)
        Phitmp=potential3D_fieldresolved(srctermall,sig0scaledall,sigPscaledall,sigHscaledall, &
                                   Vminx1,Vmaxx1,Vminx2,Vmaxx2,Vminx3,Vmaxx3, &
                                   x,flagdirich,perflag,it)
        Phiall(1:lx1,1:lx2all,1:lx3all)=Phitmp
      !else
      !  do iid=1,mpi_cfg%lid-1
      !    call mpi_send(0,1,MPI_INTEGER,iid,tag%flagdirich,MPI_COMM_WORLD)
      !  end do
      !    if (debug) print*, 'Boundary conditions too small to require solve, setting everything to zero...'
      !  Phiall=0e0_wp
      !end if
      !R------
    end if
  else   !lx2 or lx3=1 so do a field-resolved 2D solve over x1,x3
    if (debug) print *, 'Beginning field-resolved 2D solve...  Type;  ',flagdirich

    !-------
    !PRODUCE SCALED CONDUCTIVITIES TO PASS TO SOLVER, ALSO SCALED SOURCE TERM
    sig0scaled=x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h1(1:lx1,1:lx2,1:lx3)*sig0
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
    srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
    !-------

    !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
    call gather_recv(sigPscaled,tag%sigP,sigPscaledall)
    call gather_recv(sig0scaled,tag%sig0,sig0scaledall)
    call gather_recv(srcterm,tag%src,srctermall)

    !> Need to get the physical parallel conductivity so that we can convert boundary conditions for solve from current to potential
    ! Note that it is a little inefficient to have root do this calculation, but if we are in 2D it probably doesn't matter anyway...
    !print*, 'set field aligned current'
    if (flagdirich==0) then
      if (gridflag==1) then    !inverted
        Vminx1slab=-x%h1(1,1:lx2,1:lx3)*Vminx1slab/sig0(1,:,:)
        call gather_recv(Vminx1slab,tag%Vminx1,Vminx1)
      else                !non-inverted
        Vmaxx1slab=-x%h1(lx1,1:lx2,1:lx3)*Vmaxx1slab/sig0(lx1,:,:)
        call gather_recv(Vmaxx1slab,tag%Vmaxx1,Vmaxx1)
      end if
    end if

    !! EXECUTE THE SOLVE WITH MUMPS AND SCALED TERMS
    !! NOTE THE LACK OF A SPECIAL CASE HERE TO CHANGE THE POTENTIAL PROBLEM
    !! - ONLY THE HALL TERM CHANGES (SINCE RELATED TO EXB) BUT THAT DOESN'T APPEAR IN THIS EQN!
    call cpu_time(tstart)
    if (lx3all==1) then
      Phitmp=potential2D_fieldresolved(srctermall,sig0scaledall, &
                        sigPscaledall,Vminx1,Vmaxx1,Vminx2,Vmaxx2, &
                        x,flagdirich,perflag,it)
    else if (lx2all==1) then
      Phitmp=potential2D_fieldresolved(srctermall,sig0scaledall, &
                        sigPscaledall,Vminx1,Vmaxx1,Vminx3,Vmaxx3, &
                        x,flagdirich,perflag,it)
    else
      error stop '  potential_mumps.f90:  incorrect gridswap value!!!'
    end if
    Phiall(1:lx1,1:lx2all,1:lx3all)=Phitmp
    call cpu_time(tfin)
  end if
  if (debug) print *, 'MUMPS time:  ',tfin-tstart
  !!!!!!!!!

  !RADD--- ROOT NEEDS TO PUSH THE POTENTIAL BACK TO ALL WORKERS FOR FURTHER PROCESSING (BELOW)
  call bcast_send3D_ghost(Phiall,tag%Phi,Phi)

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
  !FIXME
  call pot2perpfield(Phi,x,E2,E3)

  !R-------
  !JUST TO JUDGE THE IMPACT OF MI COUPLING
  if (debug) then
    print *, 'Max integrated inertial capacitance:  ',maxval(incapintall)
    !print *, 'Max integrated Pedersen conductance (includes metric factors):  ',maxval(SigPint2all)
    !print *, 'Max integrated Hall conductance (includes metric factors):  ',minval(SigHintall), maxval(SigHintall)
    print *, 'Max E2,3 BG and response values are:  ',maxval(E02src), maxval(E03src),maxval(E2),maxval(E3)
    print *, 'Min E2,3 BG and response values are:  ',minval(E02src), minval(E03src),minval(E2),minval(E3)
    print *, 'Min/Max values of potential:  ',minval(Phi(1:lx1,1:lx2,1:lx3)),maxval(Phi(1:lx1,1:lx2,1:lx3))
    print *, 'Min/Max values of full grid potential:  ',minval(Phiall(1:lx1,1:lx2all,1:lx3all)), &
               maxval(Phiall(1:lx1,1:lx2all,1:lx3all))
  endif
  !R-------

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

  !R-------
  if (debug) then
    print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,1:lx2,1:lx3)))    !ZZZ - this rey needsz to be current at the "top"
    print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)), &
                 maxval(abs(J3pol))
    !    print *, 'Max conduction J2,3 (abs. val.) computed to be:  ',maxval(abs(J2)), &
    !                 maxval(abs(J3))
    print *, 'Max conduction J2,3  computed to be:  ',maxval(J2(1:lx1,1:lx2,1:lx3)), &
                 maxval(J3(1:lx1,1:lx2,1:lx3))
    print *, 'Min conduction J2,3  computed to be:  ',minval(J2(1:lx1,1:lx2,1:lx3)), &
                 minval(J3(1:lx1,1:lx2,1:lx3))
    print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1(1:lx1,1:lx2,1:lx3)))
  endif
  !R-------

  !-------
  !GRAND TOTAL FOR THE CURRENT DENSITY:  TOSS IN POLARIZATION CURRENT SO THAT OUTPUT FILES ARE CONSISTENT
  J1(1:lx1,1:lx2,1:lx3)=J1(1:lx1,1:lx2,1:lx3)+J1pol
  J2(1:lx1,1:lx2,1:lx3)=J2(1:lx1,1:lx2,1:lx3)+J2pol
  J3(1:lx1,1:lx2,1:lx3)=J3(1:lx1,1:lx2,1:lx3)+J3pol
  !-------
end procedure potential_root_mpi_curv

end submodule potential_root
