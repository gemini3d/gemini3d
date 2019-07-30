submodule (ionization) glow_mod

use phys_consts, only : wp

!> subroutine to allocate use-associated variables
use cglow,only: cglow_init, &
  jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst, &
  idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec, &
  iscale,jlocal,kchem,xuvfac, &
  sza,dip,efrac,ierr, &
  zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte, &
  ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion, &
  photoi,photod,phono,aglw,ecalc,zxden,zeta,zceta,eheat,vcb, &
  data_dir

implicit none

logical :: first_call = .true.

contains

module procedure glow_run

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


real, dimension(nbins) :: phitoptmp
integer :: j
character(len=1024) :: iri90_dir

data_dir    ='objects/_deps/ncarglow-src/data/'
iri90_dir   ='objects/_deps/ncarglow-src/data/iri90/'

!! Execute:  Allocate arrays in other modules (formerly in common blocks):

if(first_call) then
  first_call = .false.
  jmax=size(alt,1)
  call cglow_init
end if

!! Set electron energy grid:

call egrid (ener, del, nbins)

!! Set Maxwellian distribution into phitop array

!! Hard coded solution, future = pass ec and ef array to maxt assuming > 2 populations
phitop=0.0_wp
do j = 1, size(PhiWmWm2,1)
  !write(*,*) 'Eo and Q values: ',W0(j),PhiWmWm2(j)
  call maxt(real(PhiWmWm2(j),4),real(W0(j),4),ener,del,nbins,0,0,0,phitoptmp)
  phitop=phitop+phitoptmp
end do
!write(*,*) 'Max flux in phitop: ',maxval(pack(phitop,.true.))

!! Set variables given from GEMINI

glat = real(xlat,4)
glong = real(xlon,4)
idate = date_doy
ut = real(UTsec,4)
f107 = real(xf107,4)
f107p = real(xf107,4)
f107a = real(xf107a,4)
ap = 5.
kchem = 4.
jlocal = 0.

!! Convert densities and altitudes into

zz(:)  = real(alt(:)*1.0d2,4)
zo(:)  = real(nn(:,1)/1.0d6,4)
zo2(:) = real(nn(:,3)/1.0d6,4)
zn2(:) = real(nn(:,2)/1.0d6,4)
zno(:) = real(nn(:,6)/1.0d6,4)
zns(:) = real(nn(:,5)/1.0d6,4)
znd(:) = real(0d0,4)
ztn(:) = real(Tn(:),4)

! ZXDEN   array of excited and and/or ionized state densities at each altitude:
!           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
!           N(2D), O(1S), O(1D); cm-3
!ions (ns): 1=O+, 2=NO+, 3=N2+, 4=O2+, 5=N+, 6=H+
!neutrals (nn): O,N2,O2,H,N,NO
zxden(1, :) = 0d0
zxden(2, :) = 0d0
zxden(3, :) = real(ns(:,1)/1.0d6,4)
zxden(4, :) = real(ns(:,5)/1.0d6,4)
zxden(5, :) = real(ns(:,3)/1.0d6,4)
zxden(6, :) = real(ns(:,4)/1.0d6,4)
zxden(7, :) = real(ns(:,2)/1.0d6,4)
zxden(8, :) = 0d0
zxden(9, :) = 0d0
zxden(10,:) = 0d0
zxden(11,:) = 0d0
zxden(12,:) = 0d0
zti(:) = real((Ts(:,1)*ns(:,1)+Ts(:,2)*ns(:,2)+Ts(:,3)*ns(:,3)+Ts(:,4)*ns(:,4)+Ts(:,5)*ns(:,5)) &
        /(ns(:,1)+ns(:,2)+ns(:,3)+ns(:,4)+ns(:,5)),4)
zte(:) = real(Ts(:,7),4)
ze(:)  = real(ns(:,7)/1.0d6,4)

!! Call GLOW to calculate ionized and excited species, and airglow emission rates:

call glow

ionrate(:,1) = (real(SION(1,:),wp)+real(SION(2,:),wp)*0.3d0)*1.0d6 !O+
ionrate(:,4) = (real(SION(2,:),wp)*0.7d0)*1.0d6 !O2+
ionrate(:,3) = (real(SION(3,:),wp)*0.84d0)*1.0d6 !N2+
ionrate(:,5) = (real(SION(3,:),wp)*0.16d0)*1.0d6 !N+
ionrate(:,2) = real(0d0,wp)*1.0d6 !NO+
ionrate(:,6) = real(0d0,wp)*1.0d6 !H+
eheating = real(eheat,wp)*1.0d6
iver = real(vcb,wp)

!  do j = 1, jmax
!    IDENS(j) = sum((dflx(:,j)-uflx(:,j))*del(:))/6.241509d14
!  enddo

!  if (first_out) then
!    first_out = .false.
!    write(6,"(1x,i7,11f10.3)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec,ef1,ec1
!    write(6,"('   Z     Tn   Ti   Te      O        N2        NO      Ne(in)    Ne(out)   Ionrate     &
!          O+        O2+       NO+       N2+     EHeat    Jz')")
!    write(6,"(1x,0p,f5.1,3f6.0,1p,12e10.2)") (alt(j)/1.0d3,ztn(j),zti(j),zte(j),zo(j),zn2(j),zno(j),ze(j), &
!      ecalc(j),tir(j),PO(j),PO2(j),PNO(j),PN2(j),EHEATING(j),IDENS(j),j=1,jmax)
!  end if

end procedure glow_run

end submodule glow_mod
