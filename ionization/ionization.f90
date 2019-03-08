module ionization

use phys_consts, only: elchrg, lsp, kb, mn, re, pi, wp, lwave
use calculus, only: chapman_a
use neutral, only: Tnmsis       !we need the unperturbed msis temperatures to apply the simply chapman theory used by this module
use grid, only: curvmesh,lx1,lx2,lx3,g1,g2,g3
use timeutils, only: doy_calc
use mpi, only: MPI_COMM_WORLD,MPI_STATUS_IGNORE
use mpimod, only: myid,ierr,lid,mpi_realprec,tagTninf

implicit none
private
public :: ionrate_fang08, ionrate_glow98, eheating, photoionization

interface
module subroutine glow_run(W0,PhiWmWm2,date_doy,UTsec,xf107,xf107a,xlat,xlon,alt,nn,Tn,ns,Ts,&
  ionrate,eheating,iver)
  
real(wp), dimension(:), intent(in) :: W0,PhiWmWm2,alt,Tn
real(wp), dimension(:,:), intent(in) :: nn,ns,Ts
real(wp), dimension(:,:), intent(out) :: ionrate
real(wp), dimension(:), intent(out) :: eheating, iver
real(wp), intent(in) :: UTsec, xlat, xlon, xf107, xf107a
integer, intent(in) :: date_doy

end subroutine glow_run
end interface

! FIXME: THERE *IS* A BETTER WAY TO DO THIS
real(wp), dimension(8,4) :: Pijcoeff
data Pijcoeff(1,:) /3.49979d-1, -6.18200d-2, -4.08124d-2, 1.65414d-2/
data Pijcoeff(2,:) /5.85425d-1, -5.00793d-2, 5.69309d-2, -4.02491d-3/
data Pijcoeff(3,:) /1.69692d-1, -2.58981d-2, 1.96822d-2, 1.20505d-3/
data Pijcoeff(4,:) /-1.22271d-1, -1.15532d-2, 5.37951d-6, 1.20189d-3/
data Pijcoeff(5,:) /1.57018d0, 2.87896d-1, -4.14857d-1, 5.18158d-2/
data Pijcoeff(6,:) /8.83195d-1, 4.31402d-2, -8.33599d-2, 1.02515d-2/
data Pijcoeff(7,:) /1.90953d0, -4.74704d-2, -1.80200d-1, 2.46652d-2/
data Pijcoeff(8,:) /-1.29566d0, -2.10952d-1, 2.73106d-1,  -2.92752d-2/

contains


function photoionization(x,nn,Tn,chi,f107,f107a)

!------------------------------------------------------------
!-------COMPUTE PHOTOIONIZATION RATES PER SOLOMON ET AL, 2005
!------------------------------------------------------------

type(curvmesh), intent(in) :: x
real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: Tn
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
integer :: iid

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
Iinf=fref*(1d0+Aeuv*(0.5d0*(f107+f107a)-80d0))


!GRAVITATIONAL FIELD AND AVERAGE VALUE
g=sqrt(g1**2+g2**2+g3**2)
!    gavg=sum(g)/(lx1*lx2*lx3)    !single average value for computing colunn dens.  Interestingly this is a worker average...  Do we need root grav vars. grid mod to prevent tearing?  Should be okay as long as the grid is only sliced along the x3-dimension, BUT it isn't for simulations where arrays get permuted!!!
gavg=8d0

Tninf=maxval(Tnmsis)   !set exospheric temperature based on the max value of the background MSIS atmosphere; note this is a worker max

!both g and Tinf need to be computed as average over the entire grid...
if (myid==0) then     !root
  do iid=1,lid-1
      call mpi_recv(Tninftmp,1,mpi_realprec,iid,tagTninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      if (Tninf < Tninftmp) Tninf=Tninftmp
  end do

  do iid=1,lid-1
    call mpi_send(Tninf,1,mpi_realprec,iid,tagTninf,MPI_COMM_WORLD,ierr)
  end do  

  print *, 'Exospheric temperature used for photoionization:  ',Tninf
else                  !workders
  call mpi_send(Tninf,1,mpi_realprec,0,tagTninf,MPI_COMM_WORLD,ierr)                        !send what I think Tninf should be
  call mpi_recv(Tninf,1,mpi_realprec,0,tagTninf,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)      !receive roots decision 
end if


!O COLUMN DENSITY
H=kB*Tninf/mn(1)/gavg      !scalar scale height
bigX=(x%alt+Re)/H          !a reduced altitude
y=sqrt(bigX/2d0)*abs(cos(chi)) 
Chfn=0d0
where (chi<pi/2d0)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
!      Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*(1d0-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*(1d0+erf(y))
end where
nOcol=nn(:,:,:,1)*H*Chfn


!N2 COLUMN DENSITY
H=kB*Tninf/mn(2)/gavg     !all of these temp quantities need to be recomputed for eacb neutral species being ionized
bigX=(x%alt+Re)/H
y=sqrt(bigX/2d0)*abs(cos(chi))
Chfn=0d0
where (chi<pi/2d0)
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*(1d0+erf(y))
end where
nN2col=nn(:,:,:,2)*H*Chfn


!O2 COLUMN DENSITY
H=kB*Tninf/mn(3)/gavg
bigX=(x%alt+Re)/H
y=sqrt(bigX/2d0)*abs(cos(chi))
Chfn=0d0
where (chi<pi/2d0)    !where does work with array corresponding elements provided they are conformable
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*erfc(y)
elsewhere
  Chfn=sqrt(pi/2d0*bigX)*exp(y**2)*(1d0+erf(y))
end where
nO2col=nn(:,:,:,3)*H*Chfn


!PHOTON FLUX
do il=1,ll
  Iflux(:,:,:,il)=Iinf(il)*exp(-(sigmaO(il)*nOcol+sigmaN2(il)*nN2col+sigmaO2(il)*nO2col))
end do


!PRIMARY AND SECONDARY IONIZATION RATES
photoionization=0d0

!direct O+ production
do il=1,ll
  photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1d0+pepiO(il))
end do

!direct NO+
photoionization(:,:,:,2)=0d0

!direct N2+
do il=1,ll
  photoionization(:,:,:,3)=photoionization(:,:,:,3)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1d0+pepiN2i(il))
end do

!dissociative ionization of N2 leading to N+
do il=1,ll
  photoionization(:,:,:,5)=photoionization(:,:,:,5)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1d0+pepiN2di(il))
end do

!direct O2+
do il=1,ll
  photoionization(:,:,:,4)=photoionization(:,:,:,4)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1d0+pepiO2i(il))
end do

!dissociative ionization of O2 leading to O+
do il=1,ll
  photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1d0+pepiO2di(il))
end do

!H+ production
photoionization(:,:,:,6)=0d0


!THERE SHOULD BE SOME CODE HERE TO ZERO OUT THE BELOW-GROUND ALTITUDES.
where (photoionization<0d0)
  photoionization=0d0
end where
do isp=1,lsp-1
  phototmp=photoionization(:,:,:,isp)
  where (x%nullpts>0.9 .and. x%nullpts<1.1)
    phototmp=0d0
  end where
  photoionization(:,:,:,isp)=phototmp
end do

end function photoionization


pure function ionrate_fang08(W0,PhiWmWm2,alt,nn,Tn)

!------------------------------------------------------------
!-------COMPUTE IONIZATION RATES PER THE FANG 2008 SEMI-EMPIRICAL
!-------METHOD.
!------------------------------------------------------------

real(wp), dimension(:,:), intent(in) :: W0,PhiWmWm2

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(:,:,:), intent(in) :: alt,Tn

real(wp) :: W0keV,PhiW,alt0,deps,massden0
real(wp), dimension(1:size(nn,1)) :: H,massden,y,f,meanmass
real(wp), dimension(8) :: C

integer :: ix1,ix2,ix3,lx1,lx2,lx3,ii,ij,li,lj

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)) :: Ptot,PO,PN2,PO2

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1) :: ionrate_fang08


lx1=size(nn,1)
lx2=size(nn,2)
lx3=size(nn,3)

li=size(Pijcoeff,1)
lj=size(Pijcoeff,2)

!zero flux should really be check per field line
if ( maxval(PhiWmWm2) > 0._wp) then   !only compute rates if nonzero flux given
  !IONIZATION RATES ARE COMPUTED ON A PER-PROFILE BASIS
  deps=35d-3    !keV, kinetic energy lost per ion-electron pair produced
  do ix3=1,lx3
    do ix2=1,lx2
      !CONVERSION TO DIFFERENTIAL NUMBER FLUX
      PhiW=PhiWmWm2(ix2,ix3)*1d-3/elchrg    !from mW/m^2 to eV/m^2/s
      PhiW=PhiW/1d3/1d4    !to keV/cm^2/s
      W0keV=W0(ix2,ix3)/1d3


      !SCALE HEIGHT CALCULATION
      massden=mn(1)*nn(:,ix2,ix3,1)+mn(2)*nn(:,ix2,ix3,2)+mn(3)*nn(:,ix2,ix3,3)
      meanmass=massden/(nn(:,ix2,ix3,1)+nn(:,ix2,ix3,2)+nn(:,ix2,ix3,3))    !mean mass per particle
      H=kb*Tn(:,ix2,ix3)/meanmass/abs(g1(:,ix2,ix3))


      !Y PARAMETER
      massden=massden*1d3/1d6    !conversion of mass density from kg/m^3 to g/cm^3
      H=H*1d2                    !convert from m to cm
      y=1d0/W0keV*(massden*H/4d-6)**(0.606d0)


      !Ci COEFFS and SHAPE FUNCTION
      C=0d0
      do ii=1,li
        do ij=1,lj
          C(ii)=C(ii)+Pijcoeff(ii,ij)*log(W0keV)**(real(ij,8)-1d0)    !check whether log or log10?
        end do
        C(ii)=exp(C(ii))
      end do

      f=C(1)*y**C(2)*exp(-1*C(3)*y**C(4))+C(5)*y**C(6)*exp(-1*C(7)*y**C(8))


      !TOTAL IONIZATION RATE
      Ptot(:,ix2,ix3)=PhiW/2d0/deps/H*f*1d6   !convert to 1/m^3/s
    end do
  end do


  !NOW THAT TOTAL IONIZATION RATE HAS BEEN CALCULATED BREAK IT INTO DIFFERENT ION PRODUCTION RATES
  PO=0d0
  PN2=0d0
  PO2=0d0

  do ix3=1,lx3
    do ix2=1,lx2
      do ix1=1,lx1
        if (nn(ix1,ix2,ix3,1)+nn(ix1,ix2,ix3,2)+nn(ix1,ix2,ix3,3) > 1d-10 ) then
          PN2(ix1,ix2,ix3)=Ptot(ix1,ix2,ix3)*0.94d0*nn(ix1,ix2,ix3,2)/ &
                           (nn(ix1,ix2,ix3,3)+0.94d0*nn(ix1,ix2,ix3,2)+0.55d0*nn(ix1,ix2,ix3,1))
        end if

        if (nn(ix1,ix2,ix3,2)>1d-10) then
          PO2(ix1,ix2,ix3)=PN2(ix1,ix2,ix3)*1.07d0*nn(ix1,ix2,ix3,3)/nn(ix1,ix2,ix3,2)
          PO(ix1,ix2,ix3)=PN2(ix1,ix2,ix3)*0.59d0*nn(ix1,ix2,ix3,1)/nn(ix1,ix2,ix3,2)
        end if
      end do
    end do
  end do


  !SPLIT TOTAL IONIZATION RATE PER VALLANCE JONES, 1973
  ionrate_fang08(:,:,:,1)=PO+0.33d0*PO2
  ionrate_fang08(:,:,:,2)=0d0
  ionrate_fang08(:,:,:,3)=0.76d0*PN2
  ionrate_fang08(:,:,:,4)=0.67d0*PO2
  ionrate_fang08(:,:,:,5)=0.24d0*PN2
  ionrate_fang08(:,:,:,6)=0d0
else
  ionrate_fang08(:,:,:,:)=0d0
end if

end function ionrate_fang08


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

R=log(ns(1:lx1,1:lx2,1:lx3,lsp)/(nn(:,:,:,2)+nn(:,:,:,3)+0.1d0*nn(:,:,:,1)))
avgenergy=exp(-(12.75d0+6.941d0*R+1.166d0*R**2+0.08034d0*R**3+0.001996d0*R**4))
totionrate=sum(ionrate,4)

eheating=elchrg*avgenergy*totionrate

end function eheating


function ionrate_glow98(W0,PhiWmWm2,ymd,UTsec,f107,f107a,glat,glon,alt,nn,Tn,ns,Ts,eheating,iver)

!------------------------------------------------------------
!-------COMPUTE IONIZATION RATES USING GLOW MODEL RUN AT EACH
!-------X,Y METHOD.
!------------------------------------------------------------

real(wp), dimension(:,:,:), intent(in) :: W0,PhiWmWm2

integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec, f107, f107a
real(wp), dimension(:,:), intent(in) :: glat,glon

real(wp), dimension(:,:,:,:), intent(in) :: nn
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: ns,Ts
real(wp), dimension(:,:,:), intent(in) :: alt,Tn

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3)), intent(out) :: eheating
real(wp), dimension(1:size(nn,2),1:size(nn,3),lwave), intent(out) :: iver

real(wp), dimension(1:size(nn,1),1:size(nn,2),1:size(nn,3),lsp-1) :: ionrate_glow98

integer :: ix2,ix3,lx1,lx2,lx3,date_doy

lx1=size(nn,1)
lx2=size(nn,2)
lx3=size(nn,3)

!zero flux should really be check per field line
if ( maxval(PhiWmWm2) > 0.0_wp) then   !only compute rates if nonzero flux given
  date_doy = mod(ymd(1),100)*1000+doy_calc(ymd)

  do ix3=1,lx3
    do ix2=1,lx2
      !W0eV=W0(ix2,ix3) !Eo in eV at upper x,y locations (z,x,y) normally
     
      if ( maxval(PhiWmWm2(ix2,ix3,:)) <= 0.0_wp) then
        ionrate_glow98(:,ix2,ix3,:)=0.0_wp
        eheating(:,ix2,ix3)=0.0_wp
        iver(ix2,ix3,:)=0.0_wp
      else
        !Run GLOW here with the input parameters to obtain production rates
        !GLOW outputs ion production rates in 1/cm^3/s
        call glow_run(W0(ix2,ix3,:),PhiWmWm2(ix2,ix3,:),date_doy,UTsec,f107,f107a,glat(ix2,ix3), &
        glon(ix2,ix3),alt(:,ix2,ix3),nn(:,ix2,ix3,:),Tn(:,ix2,ix3),ns(1:lx1,ix2,ix3,:), &
        Ts(1:lx1,ix2,ix3,:),ionrate_glow98(:,ix2,ix3,:),eheating(:,ix2,ix3),iver(ix2,ix3,:))
        !write(*,*) 'glow called, max precip: ', maxval(ionrate_glow98(:,ix2,ix3,:))
      end if
    end do !Y coordinate loop
  end do !X coordinate loop
  eheating=eheating*elchrg
else
  ionrate_glow98(:,:,:,:)=0.0_wp !No Q for incoming electrons, no electron impact
  eheating(:,:,:)=0.0_wp
  iver(:,:,:)=0.0_wp
end if

end function ionrate_glow98

end module ionization
