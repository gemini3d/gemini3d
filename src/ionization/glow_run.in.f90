submodule (ionization) glow_mod

use, intrinsic :: iso_fortran_env, only: sp=>real32

!> subroutine to allocate use-associated variables
use cglow,only: cglow_init, &
  jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst, &
  idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec, &
  iscale,jlocal,kchem,xuvfac, &
  sza,dip,efrac, &
  zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte, &
  ener,del,phitop,wave1,wave2,sflux,pespec,sespec,uflx,dflx,sion, &
  photoi,photod,phono,aglw,ecalc,zxden,zeta,zceta,eheat,vcb, &
  data_dir

implicit none (type, external)

logical :: first_call = .true.

external :: egrid, maxt, glow

contains

module procedure glow_run

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon and Ben Foster, 1/15
! Stan Solomon, 12/15, 1/16
! Stan Solomon, 3/16, MPI parallel version
! Guy Grubbs II, 4/17, modified for GEMINI integration

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


real(sp), dimension(:), allocatable :: phitoptmp
integer :: j
character(len=1024) :: iri90_dir

data_dir = "@glow_data_dir@"
iri90_dir = trim(data_dir) // '/iri90/'


!! Execute:  Allocate arrays in other modules (formerly in common blocks):

if(first_call) then
  first_call = .false.
  jmax=size(alt,1)
  nbins = 190     !< MZ - note nbins no longer set to a default value in glow; user MUST assign this!!!
  call cglow_init
end if

!! Set electron energy grid:

call egrid (ener, del, nbins)

!! Set Maxwellian distribution into phitop array

!! Hard coded solution, future = pass ec and ef array to maxt assuming > 2 populations
allocate(phitoptmp(nbins))
phitoptmp=0

phitop=phitoptmp
do j = 1, size(PhiWmWm2,1)    !this index loops over population number
  call maxt(real(PhiWmWm2(j), sp),real(W0(j), sp),ener,del,nbins,0,0,0,phitoptmp)
  phitop = phitop + phitoptmp
end do


!! Set variables given from GEMINI
glat = real(xlat, sp)
glong = real(xlon, sp)
idate = date_doy
ut = real(UTsec, sp)
f107 = real(xf107, sp)
f107p = real(xf107, sp)
f107a = real(xf107a, sp)
ap = 5.
kchem = 4.
jlocal = 0

!! Convert densities and altitudes into

zz(:)  = real(alt(:)*1e2_wp, sp)
zo(:)  = real(nn(:,1)/1e6_wp, sp)
zo2(:) = real(nn(:,3)/1e6_wp, sp)
zn2(:) = real(nn(:,2)/1e6_wp, sp)
zno(:) = real(nn(:,6)/1e6_wp, sp)
zns(:) = real(nn(:,5)/1e6_wp, sp)
znd(:) = real(0, sp)
ztn(:) = real(Tn(:), sp)

! ZXDEN   array of excited and and/or ionized state densities at each altitude:
!           O+(2P), O+(2D), O+(4S), N+, N2+, O2+, NO+, N2(A), N(2P),
!           N(2D), O(1S), O(1D); cm-3
!ions (ns): 1=O+, 2=NO+, 3=N2+, 4=O2+, 5=N+, 6=H+
!neutrals (nn): O,N2,O2,H,N,NO
zxden(1, :) = 0
zxden(2, :) = 0
zxden(3, :) = real(ns(:,1)/1e6_wp, sp)
zxden(4, :) = real(ns(:,5)/1e6_wp, sp)
zxden(5, :) = real(ns(:,3)/1e6_wp, sp)
zxden(6, :) = real(ns(:, sp)/1e6_wp, sp)
zxden(7, :) = real(ns(:,2)/1e6_wp, sp)
zxden(8, :) = 0
zxden(9, :) = 0
zxden(10,:) = 0
zxden(11,:) = 0
zxden(12,:) = 0
zti(:) = real((Ts(:,1)*ns(:,1)+Ts(:,2)*ns(:,2)+Ts(:,3)*ns(:,3)+Ts(:,4)*ns(:,4)+Ts(:,5)*ns(:,5)) &
        /(ns(:,1)+ns(:,2)+ns(:,3)+ns(:,4)+ns(:,5)), sp)
zte(:) = real(Ts(:,7), sp)
ze(:)  = real(ns(:,7)/1e6_wp, sp)

!! Call GLOW to calculate ionized and excited species, and airglow emission rates:

call glow

ionrate(:,1) = (real(SION(1,:),wp)+real(SION(2,:),wp)*0.3_wp)*1e6_wp !O+
ionrate(:,4) = (real(SION(2,:),wp)*0.7_wp)*1e6_wp !O2+
ionrate(:,3) = (real(SION(3,:),wp)*0.84_wp)*1e6_wp !N2+
ionrate(:,5) = (real(SION(3,:),wp)*0.16_wp)*1e6_wp !N+
ionrate(:,2) = real(0., wp)*1e6_wp !NO+
ionrate(:,6) = real(0., wp)*1e6_wp !H+

eheating = real(eheat,wp)*1e6_wp

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

deallocate(phitoptmp)

end procedure glow_run

end submodule glow_mod
