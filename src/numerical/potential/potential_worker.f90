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

integer :: flagdirich

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: paramtrim    !to hold trimmed magnetic field

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: grad2E,grad3E    !more work arrays for pol. curr.
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: DE2Dt,DE3Dt   !pol. drift
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: J1pol,J2pol,J3pol

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E01,E02,E03   !distributed background fields
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: srcterm,divJperp
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: E1prev,E2prev,E3prev
real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: Phi

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: integrand,sigintegral    !general work array for doing integrals
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: SigPint2,SigPint3,SigHint,incapint,srctermint

real(wp), dimension(0:size(E1,1)+1,0:size(E1,2)+1,0:size(E1,3)+1) :: divtmp
!! one extra grid point on either end to facilitate derivatives
real(wp), dimension(-1:size(E1,1)+2,-1:size(E1,2)+2,-1:size(E1,3)+2) :: J1halo,J2halo,J3halo
!! haloing assumes existence of two ghost cells

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: sig0scaled,sigPscaled,sigHscaled

logical :: perflag    !MUMPS stuff

real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: Vminx1slab,Vmaxx1slab

real(wp), dimension(1:size(E1,1),1:size(E1,2),1:size(E1,3)) :: v2,v3
real(wp), dimension(1:size(E1,2),1:size(E1,3)) :: v2slab,v3slab

integer :: ix1,ix2,ix3,lx1,lx2,lx3,lx3all, ierr
integer :: idleft,idright,iddown,idup

real(wp) :: tstart,tfin

integer :: flagsolve


!SIZES - PERHAPS SHOULD BE TAKEN FROM GRID MODULE INSTEAD OF RECOMPUTED?
lx1=size(sig0,1)
lx2=size(sig0,2)
lx3=size(sig0,3)


! this should always be on by default unless the user wants to turn off and recompile; ~10% savings in mumps time *per time step*
perflag=.true.


call mpi_recv(flagdirich,1,MPI_INTEGER,0,tag%flagdirich,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
if (ierr /= 0) error stop 'dirich'

!Need to broadcast background fields from root
!Need to also broadcast x1 boundary conditions for source term calculations.
call bcast_recv(E01,tag%E01)
call bcast_recv(E02,tag%E02)
call bcast_recv(E03,tag%E03)
call bcast_recv(Vminx1slab,tag%Vminx1)
call bcast_recv(Vmaxx1slab,tag%Vmaxx1)


call potential_sourceterms(sigP,sigH,sigPgrav,sigHgrav,E02,E03,vn2,vn3,B1,x,srcterm)


!    !ZZZ - DEBUG BY GETTING THE ENTIRE SOURCETERM ARRAY
!    call gather_send(srcterm,tag%src)


!!!!!!!!
!-----AT THIS POINT WE MUST DECIDE WHETHER TO DO AN INTEGRATED SOLVE OR A 2D FIELD-RESOLVED SOLVED
!-----DECIDE BASED ON THE SIZE OF THE X2 DIMENSION
if (lx2/=1) then    !either field-resolved 3D or integrated 2D solve for 3D domain
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
    if (flagswap==1) then
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
    else
      sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
    end if
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
  if (flagswap==1) then
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)/x%h3(1:lx1,1:lx2,1:lx3)*sigP    !remember to swap 2-->3
  else
    sigPscaled=x%h1(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)/x%h2(1:lx1,1:lx2,1:lx3)*sigP
  end if
  srcterm=srcterm*x%h1(1:lx1,1:lx2,1:lx3)*x%h2(1:lx1,1:lx2,1:lx3)*x%h3(1:lx1,1:lx2,1:lx3)
  !-------

  !RADD--- NEED TO GET THE RESOLVED SOURCE TERMS AND COEFFICIENTS FROM WORKERS
  call gather_send(sigPscaled,tag%sigP)
  call gather_send(sig0scaled,tag%sig0)
  call gather_send(srcterm,tag%src)

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
!      E20all=grad3D2(-1d0*Phi0all,dx2(1:lx2))
!! causes major memory leak. maybe from arithmetic statement argument?
!! Left here as a 'lesson learned' (or is it a gfortran bug...)
!      E30all=grad3D3(-1d0*Phi0all,dx3all(1:lx3all))
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
!ADD IN BACKGROUND FIELDS BEFORE DISTRIBUTING TO WORKERS
E2=E2+E02
E3=E3+E03
!--------


!--------
!COMPUTE TIME DERIVATIVE NEEDED FOR POLARIZATION CURRENT.  ONLY DO THIS IF WE HAVE SPECIFIC NONZERO INERTIAL CAPACITANCE
!if (maxval(incap) > 0._wp) then
if (cfg%flagcap/=0) then
  !differentiate E2 in x2
  J1halo(1:lx1,1:lx2,1:lx3)=E2
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D2(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad2E=divtmp(1:lx1,1:lx2,1:lx3)

  !differentiate E2 in x3
  J1halo(1:lx1,1:lx2,1:lx3)=E2
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)    !likely doesn't need to be haloed again
  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  !compute total derivative in x2
  DE2Dt=(E2-E2prev)/dt+v2*grad2E+v3*grad3E

  !differentiate E3 in x2
  J1halo(1:lx1,1:lx2,1:lx3)=E3
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)
  divtmp=grad3D2(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad2E=divtmp(1:lx1,1:lx2,1:lx3)

  !differentiate E3 in x3
  J1halo(1:lx1,1:lx2,1:lx3)=E3
  call halo_pot(J1halo,tag%J1,x%flagper,.false.)    !maybe don't need to halo again???
  divtmp=grad3D3(J1halo(0:lx1+1,0:lx2+1,0:lx3+1),x,0,lx1+1,0,lx2+1,0,lx3+1)
  grad3E=divtmp(1:lx1,1:lx2,1:lx3)

  !x3 total derivative
  DE3Dt=(E3-E3prev)/dt+v2*grad2E+v3*grad3E

  !convert derivative into polarization current density
  J1pol=0d0
  J2pol=incap*DE2Dt
  J3pol=incap*DE3Dt
else       !pure electrostatic solve was done
  DE2Dt=0d0
  DE3Dt=0d0
  J1pol=0d0
  J2pol=0d0
  J3pol=0d0
end if
!--------

!-------
if (flagswap==1) then
  J2=sigP*E2+sigH*E3    !BG field already added to E above
  J3=-1*sigH*E2+sigP*E3
else
  J2=sigP*E2-sigH*E3    !BG field already added to E above
  J3=sigH*E2+sigP*E3
end if

!WHAT I THINK THE NEUTRAL WIND CURRENTS SHOULD BE IN 2D
if (flagswap==1) then
  J2=J2-sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)+sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
else
  J2=J2+sigP*vn3*B1(1:lx1,1:lx2,1:lx3)+sigH*vn2*B1(1:lx1,1:lx2,1:lx3)
  J3=J3+sigH*vn3*B1(1:lx1,1:lx2,1:lx3)-sigP*vn2*B1(1:lx1,1:lx2,1:lx3)
end if
!-------


!!!!!!!!
!NOW DEAL WITH THE PARALLEL FIELDS AND ALL CURRENT
if (lx2/=1 .and. cfg%potsolve==1) then    !we did a field-integrated solve above
  !-------
  !NOTE THAT A DIRECT E1ALL CALCULATION WILL GIVE ZERO, SO USE INDIRECT METHOD, AS FOLLOWS
  J1=0d0

  if (cfg%flagJpar) then    !user can elect to not compute parallel current; in some low-res cases is it too prone to artifacts to use reliably
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
      !! - they have to  be loaded into VVmaxx1 and Vminx1.
      !! For numerical purposes we prefer to integrate from the location of nonzero current (usually highest altitude in open grid).
      if (gridflag==0) then     !closed dipole grid, really would be best off integrating from the source hemisphere
  !          if (debug) print *,  'Closed dipole grid; integration starting at max x1...', minval(Vmaxx1slab), &
  !                         maxval(Vmaxx1slab)
        if (cfg%sourcemlat>=0d0) then    !integrate from northern hemisphere
  !            if (debug) print *, 'Source is in northern hemisphere (or there is no source)...'
          J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
          do ix1=1,lx1
            J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                             (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
          end do
        else
  !            if (debug) print *, 'Source in southern hemisphere...'
          J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
          do ix1=1,lx1
            J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                             (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
          end do
        end if
      elseif (gridflag==1) then    !this would be an inverted grid, this max altitude corresponds to the min value of x1
  !          if (debug) print *,  'Inverted grid; integration starting at min x1...',minval(Vminx1slab), maxval(Vminx1slab)
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vminx1slab-J1(ix1,:,:))
        end do
      else        !minx1 is at teh bottom of the grid to integrate from max x1
  !          if (debug) print *,  'Non-inverted grid; integration starting at max x1...', minval(Vmaxx1slab),  maxval(Vmaxx1slab)
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current, go from maxval(x1) to location of interest
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (x%h2(1,1:lx2,1:lx3)*x%h3(1,1:lx2,1:lx3)*Vmaxx1slab+J1(ix1,:,:))
        end do
      end if
  !        if (gridflag==2) then
           !! for a cartesian grid in the northern hemisphere (assumed) we have the x1-direction being against the magnetic field...
  !          J1=-1d0*J1     !ZZZ - very questionable
  !        end if
    else
      !! Dirichlet conditions - we need to integrate from the ***lowest altitude***
      !! (where FAC is known to be zero, note this is not necessarilty the logical bottom of the grid), upwards (to where it isn't)
      if (gridflag/=2) then    !inverted grid (logical top is the lowest altitude)
  !          if (debug) print *, 'Inverted grid detected - integrating logical top downward to compute FAC...'
        J1=integral3D1_curv_alt(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (J1(ix1,:,:))    !FAC AT TOP ASSUMED TO BE ZERO
        end do
      else      !non-inverted grid (logical bottom is the lowest altitude - so integrate normy)
  !          if (debug) print *, 'Non-inverted grid detected - integrating logical bottom to top to compute FAC...'
        J1=integral3D1(divJperp,x,1,lx1)    !int divperp of BG current
        do ix1=1,lx1
          J1(ix1,:,:)=1d0/x%h2(ix1,1:lx2,1:lx3)/x%h3(ix1,1:lx2,1:lx3)* &
                           (-1d0*J1(ix1,:,:))    !FAC AT THE BOTTOM ASSUMED TO BE ZERO
        end do
      end if
    end if
    E1=J1/sig0
    !-------
  end if !flagJpar
else   !we resolved the field line (either 2D solve or full 3D) so just differentiate normally

  !-------
  Phi=-1d0*Phi
  E1=grad3D1(Phi,x,1,lx1,1,lx2,1,lx3)    !no haloing required since x1-derivative
  Phi=-1d0*Phi
  J1=sig0*E1
  !-------

end if
!!!!!!!!!


!    !R-------
if (debug) then
!    print *, 'Max topside FAC (abs. val.) computed to be:  ',maxval(abs(J1(1,:,:)))
!! ZZZ - this rey needsz to be current at the "top"
!    print *, 'Max polarization J2,3 (abs. val.) computed to be:  ',maxval(abs(J2pol)),  maxval(abs(J3pol))
!    print *, 'Max conduction J2,3  computed to be:  ',maxval(J2), maxval(J3)
!    print *, 'Min conduction J2,3  computed to be:  ',minval(J2),  minval(J3)
!    print *, 'Max conduction J1 (abs. val.) computed to be:  ',maxval(abs(J1))
!    print *, 'flagswap:  ',flagswap
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
