subroutine glow_run(Eo,Q,alt,Nn,Tn,Ns,Ts,PO,PNO,PN2,PO2,PH,EHEATING,IDENS)

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon and Ben Foster, 1/15
! Stan Solomon, 12/15, 1/16
! Stan Solomon, 3/16, MPI parallel version
! Guy Grubbs II, 4/17, modified for GEMINI integration

! Main multi-processor driver for the GLOW model.
! Uses TIE-GCM history files or MSIS/IRI for input.
! Requires MPI and netCDF libraries

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.
! For definitions of TGCM input variables see module READTGCM
! For definitions of output arrays see module OUTPUT

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major` species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission


  use cglow,only: cglow_init,cglow_reset      ! subroutine to allocate use-associated variables
  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst
  use cglow,only: idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,ef1,ec1
  use cglow,only: iscale,jlocal,kchem,xuvfac
  use cglow,only: sza,dip,efrac,ierr
  use cglow,only: zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte
  use cglow,only: ener,del,phitop,phitop2,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion
  use cglow,only: photoi,photod,phono,aglw,tei,tpi,tir,ecalc,zxden,zeta,zceta,eheat,data_dir
  use cglow,only: first_run

  implicit none
  
  character(len=1024) :: iri90_dir

  real(8), dimension(:), intent(in) :: Eo,Q 
  real(8), dimension(:), intent(in) :: alt,Tn           ! glow height coordinate in km (jmax)
  real(8), dimension(:,:), intent(in) :: Nn
  real(8), dimension(:,:), intent(in) :: Ns,Ts
  real(8), dimension(:), intent(out) :: PO,PN2,PO2,PNO,PH,EHEATING,IDENS
  
  integer :: j = 0
  logical :: first_out = .true.
!  namelist /glow_input/ &
!    indate,utstart,utstep,utstop,nlat_msis,nlon_msis,f107a,f107,f107p,ap, &
!    iscale,jlocal,kchem,xuvfac,ef,ec,itail,fmono,emono, &
!    tgcm_ncfile,iri90_dir,jmax,glow_ncfile, &
!    start_mtime,stop_mtime,data_dir,writelbh,writered
  data_dir    ='data/'
  iri90_dir   ='data/iri90'
! Execute:
! Allocate arrays in other modules (formerly in common blocks):
!
  if(first_run) then
    first_run = .false.
    jmax=size(alt,1)
    call cglow_init
    call cglow_reset
  else
    call cglow_reset
  end if
!
! Set electron energy grid:
!
  call egrid (ener, del, nbins)
!
! Set Maxwellian distribution into phitop array
!
! Hard coded solution, future = pass ec and ef array to maxt assuming > 2
! populations
  ec=Eo(1)
  ef=Q(1)
  ec1=Eo(2)
  ef1=Q(2)
  call maxt2 (ef,ec,ef1,ec1,ener,del,nbins,phitop,phitop2)
  phitop=phitop+phitop2
!
! Set variables for Rocket B launch and given densities from GEMINI
!
  glat = 67.01
  glong = 213.58
  idate = 2017061
  ut = 28800.
  f107 = 79.
  f107p = 79.
  f107a = 79.
  ap = 5.
  kchem = 4.
  jlocal = 0.
!
! Convert densities and altitudes into 
!
  zz(:)  = alt(:)*1.0d2
  zo(:)  = NN(:,1)/1.0d6
  zo2(:) = NN(:,3)/1.0d6
  zn2(:) = NN(:,2)/1.0d6
  zno(:) = NN(:,6)/1.0d6
  zns(:) = NN(:,5)/1.0d6
  znd(:) = 0d0 
  ztn(:) = Tn(:)

! ZXDEN   array of excited and and/or ionized state densities at each altitude:
!           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
!           N(2D), O(1S), O(1D); cm-3
!ions (ns): 1=O+, 2=NO+, 3=N2+, 4=O2+, 5=N+, 6=H+
  zxden(1, :) = 0d0
  zxden(2, :) = 0d0
  zxden(3, :) = Ns(:,1)/1.0d6
  zxden(4, :) = Ns(:,5)/1.0d6
  zxden(5, :) = Ns(:,3)/1.0d6
  zxden(6, :) = Ns(:,4)/1.0d6
  zxden(7, :) = Ns(:,2)/1.0d6
  zxden(8, :) = 0d0
  zxden(9, :) = 0d0 
  zxden(10,:) = 0d0
  zxden(11,:) = 0d0
  zxden(12,:) = 0d0
  zti(:) = (Ts(:,1)*Ns(:,1)+Ts(:,2)*Ns(:,2)+Ts(:,3)*Ns(:,3)+Ts(:,4)*Ns(:,4)+Ts(:,5)*Ns(:,5)) &
          /(Ns(:,1)+Ns(:,2)+Ns(:,3)+Ns(:,4)+Ns(:,5))
  zte(:) = Ts(:,7)
  ze(:)  = Ns(:,7)/1.0d6
!
! Call GLOW to calculate ionized and excited species, and airglow emission rates:
!
  call glow
  
  PO(:) = sion(1,:)*1.0d6
  PO2(:) = sion(2,:)*1.0d6
  PN2(:) = sion(3,:)*1.0d6
  PNO(:) = 0d0 !Possibly edit later
  PH(:) = 0d0 !Possibly edit later
  EHEATING(:) = eheat(:)*1.0d6
  
  do j = 1, jmax
    IDENS(j) = sum((dflx(:,j)-uflx(:,j))*del(:))/6.241509d14
  enddo
  
!  if (first_out) then
!    first_out = .false.
!    write(6,"(1x,i7,11f10.3)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,ef1,ec1
!    write(6,"('   Z     Tn   Ti   Te      O        N2        NO      Ne(in)    Ne(out)   Ionrate     &
!          O+        O2+       NO+       N2+     EHeat    Jz')")
!    write(6,"(1x,0p,f5.1,3f6.0,1p,12e10.2)") (alt(j)/1.0d3,ztn(j),zti(j),zte(j),zo(j),zn2(j),zno(j),ze(j), &
!      ecalc(j),tir(j),PO(j),PO2(j),PNO(j),PN2(j),EHEATING(j),IDENS(j),j=1,jmax)
!  end if

end subroutine glow_run
