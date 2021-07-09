module ionization

use phys_consts, only: elchrg, lsp, kb, mn, re, pi, wp, lwave, debug
use neutral, only: Tnmsis
use ionize_fang, only: fang2008, fang2010
!! we need the unperturbed msis temperatures to apply the simple chapman theory used by this module
use grid, only: lx1,lx2,lx3,g1,g2,g3
use meshobj, only: curvmesh
use timeutils, only: ymd2doy
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi, MPI_COMM_WORLD,MPI_STATUS_IGNORE

implicit none (type, external)
private
public :: ionrate_fang, ionrate_glow98, eheating, photoionization

external :: mpi_send, mpi_recv

interface
module subroutine glow_run(W0,PhiWmWm2,date_doy,UTsec,xf107,xf107a,xlat,xlon,alt,nn,Tn,ns,Ts,&
  ionrate,eheating,iver)

real(wp), dimension(:), intent(in) :: W0,PhiWmWm2,alt,Tn
real(wp), dimension(:,:), intent(in) :: nn,ns,Ts
real(wp), dimension(:,:), intent(inout) :: ionrate
!! intent(out)
real(wp), dimension(:), intent(inout) :: eheating, iver
!! intent(out)
real(wp), intent(in) :: UTsec, xlat, xlon, xf107, xf107a
integer, intent(in) :: date_doy

end subroutine glow_run
end interface

contains


function photoionization(x,nn,chi,f107,f107a)

!------------------------------------------------------------
!-------COMPUTE PHOTOIONIZATION RATES PER SOLOMON ET AL, 2005
!------------------------------------------------------------

class(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:,:), intent(in) :: nn
!real(wp), dimension(:,:,:), intent(in) :: Tn
real(wp), dimension(:,:,:), intent(in) :: chi
real(wp), intent(in) :: f107,f107a

integer, parameter :: ll=22     !number of wavelength bins
integer :: il,isp
real(wp), dimension(ll) :: lambda1,lambda2,fref,Aeuv,sigmaO,sigmaN2,sigmaO2
real(wp), dimension(ll) :: brN2i,brN2di,brO2i,brO2di,pepiO,pepiN2i,pepiN2di,pepiO2i,pepiO2di
real(wp), dimension(ll) :: Iinf
real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: g,bigX,y,Chfn
real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol,nN2col,nO2col
real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: phototmp
real(wp) :: gavg,H,Tninf
real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),ll) :: Iflux

real(wp) :: Tninftmp
integer :: iid, ierr

real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),lsp-1) :: photoionization    !don't need a separate rate for electrons


!WAVELENGTH BIN BEGINNING AND END (THIS IDEALLY WOULD BE DATA STATEMENTS OR SOME KIND OF STRUCTURE THAT DOESN'T GET REASSIGNED AT EVERY CALL).  Actually all of these array assignments are static...
lambda1=[0.05, 0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 65.0, &
    79.8, 79.8, 79.8, 91.3, 91.3, 91.3, 97.5, 98.7, 102.7]*1e-9
lambda2=[0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 79.8, 79.8, &
     91.3, 91.3, 91.3, 97.5, 97.5, 97.5, 98.7, 102.7, 105.0]*1e-9


!EUVAC FLUX VALUES
fref=[5.01e1, 1e4, 2e6, 2.85e7, 5.326e8, 1.27e9, 5.612e9, 4.342e9, 8.380e9, &
    2.861e9, 4.83e9, 1.459e9, 1.142e9, 2.364e9, 3.655e9, 8.448e8, 3.818e8, &
    1.028e9, 7.156e8, 4.482e9, 4.419e9, 4.235e9]*1e4        !convert to m^-2 s^-1
Aeuv=[6.24e-1, 3.71e-1, 2e-1, 6.247e-2, 1.343e-2, 9.182e-3, 1.433e-2, 2.575e-2, &
    7.059e-3, 1.458e-2, 5.857e-3, 5.719e-3, 3.680e-3, 5.310e-3, 5.261e-3, 5.437e-3, &
    4.915e-3, 4.995e-3, 4.422e-3, 3.950e-3, 5.021e-3, 4.825e-3]


!TOTAL ABSORPTION CROSS SECTIONS
sigmaO=[0.0023, 0.0170, 0.1125, 0.1050, 0.3247, 1.3190, 3.7832, 6.0239, &
     7.7205, 10.7175, 13.1253, 8.5159, 4.7889, 3.0031, 4.1048, 3.7947, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0]*1e-18*1e-4         !convert to m^2
sigmaN2=[0.0025, 0.0201, 0.1409, 1.1370, 0.3459, 1.5273, 5.0859, 9.9375, &
    11.7383, 19.6514, 23.0931, 23.0346, 54.5252, 2.1434, 13.1062, 71.6931, &
    2.1775, 14.4390, 115.257, 2.5465, 0.0, 0.0]*1e-18*1e-4
sigmaO2=[0.0045, 0.034, 0.2251, 0.2101, 0.646, 2.6319, 7.6283, 13.2125, &
    16.8233, 20.3066, 27.0314, 23.5669, 24.9102, 10.4980, 10.9075, 13.3122, &
    13.3950, 14.4042, 32.5038, 18.7145, 1.6320, 1.15]*1e-18*1e-4


!BRANCHING RATIOS
brN2i=[0.040,0.040,0.040,0.040, 0.717, 0.751, 0.747, 0.754, 0.908, 0.996, 1.0, 0.679,  &
    0.429, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
brN2di=[0.96, 0.96,0.96,0.96,0.282, 0.249, 0.253, 0.246, 0.093, 0.005, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
brO2i=[0.0, 0.0, 0.0, 0.0, 0.108, 0.347, 0.553, 0.624, 0.649, 0.759, 0.874, 0.672, 0.477, &
    0.549, 0.574, 0.534, 0.756, 0.786, 0.620, 0.830, 0.613, 0.0]
brO2di=[1.0, 1.0, 1.0, 1.0, 0.892, 0.653, 0.447, 0.376, 0.351, 0.240, 0.108, 0.001, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


!PHOTOELECTRON TO DIRECT PRODUCTION RATIOS
pepiO=[217.12, 50.593, 23.562, 71.378, 4.995, 2.192, 1.092, 0.694, 0.418, &
    0.127, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pepiN2i=[263.99, 62.57, 25.213, 8.54, 6.142, 2.288, 0.786, 0.324, 0.169, 0.031, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pepiN2di=[78.674, 18.310, 6.948, 2.295, 1.647, 0.571, 0.146, 0.037, 0.008, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pepiO2i=[134.69, 32.212, 13.309, 39.615, 2.834, 1.092, 0.416, 0.189, 0.090, 0.023, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
pepiO2di=[76.136, 17.944, 6.981, 20.338, 1.437, 0.521, 0.163, 0.052, 0.014, 0.001, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


!IRRADIANCE ACCORDING TO [RICHARDS, 1994]
Iinf=fref*(1 + Aeuv*(0.5_wp*(f107+f107a)-80._wp))


!GRAVITATIONAL FIELD AND AVERAGE VALUE
g=sqrt(g1**2+g2**2+g3**2)
!    gavg=sum(g)/(lx1*lx2*lx3)    !single average value for computing column dens.  Interestingly this is a worker average...  Do we need root grav vars. grid mod to prevent tearing?  Should be okay as long as the grid is only sliced along the x3-dimension, BUT it isn't for simulations where arrays get permuted!!!
gavg=8._wp

Tninf=maxval(Tnmsis)   !set exospheric temperature based on the max value of the background MSIS atmosphere; note this is a worker max

!both g and Tinf need to be computed as average over the entire grid...
if (mpi_cfg%myid==0) then     !root
  ierr=0
  do iid=1,mpi_cfg%lid-1
      call mpi_recv(Tninftmp,1,mpi_realprec,iid,tag%Tninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      if (Tninf < Tninftmp) Tninf=Tninftmp
  end do
  if (ierr /= 0) error stop 'root failed to mpi_recv Tninf'

  ierr=0
  do iid=1,mpi_cfg%lid-1
    call mpi_send(Tninf,1,mpi_realprec,iid,tag%Tninf,MPI_COMM_WORLD,ierr)
  end do
  if (ierr /= 0) error stop 'root failed to mpi_send Tninf'

  if (debug) print *, 'Exospheric temperature used for photoionization:  ',Tninf
else                  !workders
  call mpi_send(Tninf,1,mpi_realprec,0,tag%Tninf,MPI_COMM_WORLD,ierr)                        !send what I think Tninf should be
  if (ierr /= 0) error stop 'worker failed to mpi_send Tninf'
  call mpi_recv(Tninf,1,mpi_realprec,0,tag%Tninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)      !receive roots decision
  if (ierr /= 0) error stop 'worker failed to mpi_recv Tninf'
end if



!O COLUMN DENSITY
H=kB*Tninf/mn(1)/gavg      !scalar scale height
bigX=(x%alt+Re)/H          !a reduced altitude
y=sqrt(bigX/2._wp)*abs(cos(chi))
Chfn=0
where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
!      Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1._wp-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
end where
nOcol=nn(:,:,:,1)*H*Chfn


!N2 COLUMN DENSITY
H=kB*Tninf/mn(2)/gavg     !all of these temp quantities need to be recomputed for eacb neutral species being ionized
bigX=(x%alt+Re)/H
y=sqrt(bigX/2._wp)*abs(cos(chi))
Chfn=0
where (chi<pi/2._wp)
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
end where
nN2col=nn(:,:,:,2)*H*Chfn


!O2 COLUMN DENSITY
H=kB*Tninf/mn(3)/gavg
bigX=(x%alt+Re)/H
y=sqrt(bigX/2._wp)*abs(cos(chi))
Chfn=0
where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
end where
nO2col=nn(:,:,:,3)*H*Chfn


!PHOTON FLUX
do il=1,ll
  Iflux(:,:,:,il)=Iinf(il)*exp(-(sigmaO(il)*nOcol+sigmaN2(il)*nN2col+sigmaO2(il)*nO2col))
end do


!PRIMARY AND SECONDARY IONIZATION RATES
photoionization=0

!direct O+ production
do il=1,ll
  photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
end do

!direct NO+
photoionization(:,:,:,2) = 0

!direct N2+
do il=1,ll
  photoionization(:,:,:,3)=photoionization(:,:,:,3)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
end do

!dissociative ionization of N2 leading to N+
do il=1,ll
  photoionization(:,:,:,5)=photoionization(:,:,:,5)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
end do

!direct O2+
do il=1,ll
  photoionization(:,:,:,4)=photoionization(:,:,:,4)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
end do

!dissociative ionization of O2 leading to O+
do il=1,ll
  photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
end do

!H+ production
photoionization(:,:,:,6) = 0


!THERE SHOULD BE SOME CODE HERE TO ZERO OUT THE BELOW-GROUND ALTITUDES.
where (photoionization < 0)
  photoionization = 0
end where
do isp=1,lsp-1
  phototmp=photoionization(:,:,:,isp)
!  where (x%nullpts>0.9 .and. x%nullpts<1.1)
  where(x%nullpts)
    phototmp=0
  end where
  photoionization(:,:,:,isp) = phototmp
end do

end function photoionization


pure function ionrate_fang(W0, PhiWmWm2, alt, nn, Tn, flag_fang)

real(wp), dimension(:,:), intent(in) :: W0,PhiWmWm2

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: alt,Tn
integer, intent(in) :: flag_fang

real(wp) :: W0keV,PhiW
real(wp), dimension(1:size(nn,1)) :: massden,meanmass

integer :: ix2,ix3,lx2,lx3

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: Ptot,PO,PN2,PO2

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1) :: ionrate_fang


lx2=size(nn,2)
lx3=size(nn,3)

!IONIZATION RATES ARE COMPUTED ON A PER-PROFILE BASIS

!zero flux should really be check per field line
if ( maxval(PhiWmWm2) > 0) then   !only compute rates if nonzero flux given

  do ix3=1,lx3
    do ix2=1,lx2
      !CONVERSION TO DIFFERENTIAL NUMBER FLUX
      PhiW=PhiWmWm2(ix2,ix3)*1e-3_wp/elchrg    !from mW/m^2 to eV/m^2/s
      PhiW=PhiW/1e3_wp/1e4_wp    !to keV/cm^2/s
      W0keV=W0(ix2,ix3)/1e3_wp

      massden=mn(1)*nn(:,ix2,ix3,1)+mn(2)*nn(:,ix2,ix3,2)+mn(3)*nn(:,ix2,ix3,3)
      !! mass densities are [kg m^-3] as per neutral/neutral.f90 "call meters(.true.)" for MSIS.
      meanmass=massden/(nn(:,ix2,ix3,1)+nn(:,ix2,ix3,2)+nn(:,ix2,ix3,3))
      !! mean mass per particle [kg]

      !> TOTAL IONIZATION RATE
      !! [cm^-3 s^-1] => [m^-3 s^-1]
      select case (flag_fang)
      case (8, 2008)
        Ptot(:,ix2,ix3) = fang2008(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3)) * 1e6_wp
      case (10, 2010)
        Ptot(:,ix2,ix3) = fang2010(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3)) * 1e6_wp
      case default
        error stop 'ERROR:ionization:ionrate_fang: unknown flag_fang'
      end select
    end do
  end do


  !NOW THAT TOTAL IONIZATION RATE HAS BEEN CALCULATED BREAK IT INTO DIFFERENT ION PRODUCTION RATES
  PO = 0
  PN2 = 0
  PO2 = 0

  where (nn(:,:,:,1) + nn(:,:,:,2) + nn(:,:,:,3) > 1e-10_wp )
          PN2 = Ptot * 0.94_wp * nn(:,:,:,2) / &
                           (nn(:,:,:,3) + 0.94_wp*nn(:,:,:,2) + 0.55_wp * nn(:,:,:,1))

  endwhere

  where (nn(:,:,:,2) > 1e-10_wp)
    PO2 = PN2 * 1.07_wp * nn(:,:,:,3) / nn(:,:,:,2)
    PO = PN2 * 0.59_wp * nn(:,:,:,1) / nn(:,:,:,2)
  endwhere



  !SPLIT TOTAL IONIZATION RATE PER VALLANCE JONES, 1973
  ionrate_fang(:,:,:,1) = PO + 0.33_wp * PO2
  ionrate_fang(:,:,:,2) = 0
  ionrate_fang(:,:,:,3) = 0.76_wp * PN2
  ionrate_fang(:,:,:,4) = 0.67_wp * PO2
  ionrate_fang(:,:,:,5) = 0.24_wp * PN2
  ionrate_fang(:,:,:,6) = 0
else
  ionrate_fang(:,:,:,:) = 0
end if

end function ionrate_fang


pure function eheating(nn,Tn,ionrate,ns)

!------------------------------------------------------------
!-------COMPUTE SWARTZ AND NISBET, (1973) ELECTRON HEATING RATES.
!-------ION ARRAYS (EXCEPT FOR RATES) ARE EXPECTED TO INCLUDE
!-------GHOST CELLS.
!------------------------------------------------------------

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn
real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1), intent(in) :: ionrate
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns    !includes ghost cells

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: totionrate,R,avgenergy
integer :: lx1,lx2,lx3

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: eheating

lx1=size(nn,1)
lx2=size(nn,2)
lx3=size(nn,3)

R=log(ns(1:lx1,1:lx2,1:lx3,lsp)/(nn(:,:,:,2)+nn(:,:,:,3)+0.1_wp*nn(:,:,:,1)))
avgenergy=exp(-(12.75_wp+6.941_wp*R+1.166_wp*R**2+0.08034_wp*R**3+0.001996_wp*R**4))
totionrate=sum(ionrate,4)

eheating=elchrg*avgenergy*totionrate

end function eheating


subroutine ionrate_glow98(W0,PhiWmWm2,ymd,UTsec,f107,f107a,glat,glon,alt,nn,Tn,ns,Ts, &
                               eheating, iver, ionrate)

!! COMPUTE IONIZATION RATES USING GLOW MODEL RUN AT EACH
!! X,Y METHOD.

real(wp), dimension(:,:,:), intent(in) :: W0,PhiWmWm2

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec, f107, f107a
real(wp), dimension(:,:), intent(in) :: glat,glon

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
real(wp), dimension(:,:,:), intent(in) :: alt,Tn

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)), intent(inout) :: eheating
!! intent(out)
real(wp), dimension(1:size(nn,2),1:size(nn,3),lwave), intent(inout) :: iver
!! intent(out)
real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1), intent(inout) :: ionrate
!! intent(out)

integer :: ix2,ix3,lx1,lx2,lx3,date_doy

lx1=size(nn,1)
lx2=size(nn,2)
lx3=size(nn,3)

!! zero flux should really be checked per field line
if ( maxval(PhiWmWm2) > 0) then   !only compute rates if nonzero flux given

  date_doy = modulo(ymd(1), 100)*1000 + ymd2doy(ymd(1), ymd(2), ymd(3))
  !! date in format needed by GLOW (yyddd)
  do ix3=1,lx3
    do ix2=1,lx2
      !W0eV=W0(ix2,ix3) !Eo in eV at upper x,y locations (z,x,y) normally

      if ( maxval(PhiWmWm2(ix2,ix3,:)) <= 0) then    !only compute rates if nonzero flux given *here* (i.e. at this location)
        ionrate(:,ix2,ix3,:) = 0
        eheating(:,ix2,ix3) = 0
        iver(ix2,ix3,:) = 0
      else
        !Run GLOW here with the input parameters to obtain production rates
        !GLOW outputs ion production rates in [cm^-3 s^-1]
        call glow_run(W0(ix2,ix3,:), PhiWmWm2(ix2,ix3,:), &
          date_doy, UTsec, f107, f107a, glat(ix2,ix3), glon(ix2,ix3), alt(:,ix2,ix3), &
          nn(:,ix2,ix3,:),Tn(:,ix2,ix3), ns(1:lx1,ix2,ix3,:), Ts(1:lx1,ix2,ix3,:), &
          ionrate(:,ix2,ix3,:), eheating(:,ix2,ix3), iver(ix2,ix3,:))
!        print*, 'glow called, max ionization rate: ', maxval(ionrate(:,ix2,ix3,:))
!        print*, 'max iver:  ',maxval(iver(ix2,ix3,:))
!        print*, 'max W0 and Phi:  ',maxval(W0(ix2,ix3,:)),maxval(PhiWmWm2(ix2,ix3,:))
      end if
    end do !Y coordinate loop
  end do !X coordinate loop
  eheating=eheating*elchrg
else
  ionrate(:,:,:,:)=0 !No Q for incoming electrons, no electron impact
  eheating(:,:,:)=0
  iver(:,:,:)=0
end if

end subroutine ionrate_glow98

end module ionization
