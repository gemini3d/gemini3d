module ionization

use gemini3d_config, only: gemini_cfg
use phys_consts, only: elchrg, lsp, kb, mn, re, pi, wp, lwave, debug
use ionize_fang, only: fang2008, fang2010, fang2010_spectrum
!! we need the unperturbed msis temperatures to apply the simple chapman theory used by this module
use grid, only: lx1,lx2,lx3
use meshobj, only: curvmesh
use timeutils, only: ymd2doy

implicit none (type, external)
private
public :: ionrate_fang, ionrate_glow98, eheating, photoionization

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
  function photoionization(cfg,t,ymd,UTsec,x,nn,chi,f107,f107a,gavg,Tninf,Iinf)
    !------------------------------------------------------------
    !-------COMPUTE PHOTOIONIZATION RATES PER SOLOMON ET AL, 2005
    !------------------------------------------------------------
    type(gemini_cfg), intent(in) :: cfg
    real(wp), intent(in) :: t
    integer, intent(in), dimension(3) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(in) :: nn
    !real(wp), dimension(:,:,:), intent(in) :: Tn
    real(wp), dimension(:,:,:), intent(in) :: chi
    real(wp), parameter :: chi0 = pi/2._wp
    real(wp), parameter :: dchi = 5._wp * pi / 180._wp   ! 5° smooth transition
    real(wp) :: w
    real(wp) :: Fchap, Fnight
    real(wp), intent(in) :: f107,f107a
    real(wp), intent(in) :: gavg,Tninf
    real(wp), dimension(:,:,:,:), intent(in) :: Iinf
    integer, parameter :: ll=22     !number of wavelength bins (daytime bins only)
    integer, parameter :: llnight = 4 !(nighttime bins only)
    integer :: il,isp,ix1,ix2,ix3
    ! Daytime EUVAC bins remain 22 as Iinf and file-based solar-flux inputs use 22 bins.
    ! Nighttime lines are handled independently below.
    character(len=8), parameter :: night_line(llnight) = &
      [character(len=8) :: "heii", "hei", "lybeta", "lyalpha"]
    real(wp), dimension(ll) :: lambda1,lambda2,sigmaO,sigmaN2,sigmaO2
    real(wp), dimension(llnight) :: lambda1_night,lambda2_night
    real(wp), dimension(llnight) :: sigiO_night,sigiN2_night,sigiO2_night
    real(wp), dimension(llnight) :: sigaO_night,sigaN2_night,sigaO2_night
    real(wp), dimension(llnight) :: brN2i_night,brN2di_night,brO2i_night,brO2di_night
    !From Strobel 1980 & Kirby 1979 --> ionization and absorption cross sections for night lines 
    !(convert 10^-18 cm^2 to m^2) 
    real(wp), parameter :: cs = 1e-22_wp
    
    !Ionization cross sections (for production)
    ! He II (304 Å)
    real(wp), parameter :: sigiO_heii  = 9.7_wp  * cs
    real(wp), parameter :: sigiN2_heii = 11.6_wp * cs
    real(wp), parameter :: sigiO2_heii = 16.0_wp * cs

    ! He I (584 Å)
    real(wp), parameter :: sigiO_hei  = 12.2_wp * cs
    real(wp), parameter :: sigiN2_hei = 23.2_wp * cs
    real(wp), parameter :: sigiO2_hei = 22.0_wp * cs

    ! Ly-beta (1026 Å)
    real(wp), parameter :: sigiO_lyb  = 0._wp * cs
    real(wp), parameter :: sigiN2_lyb = 0._wp * cs
    real(wp), parameter :: sigiO2_lyb = 1.0_wp * cs

    ! Ly-alpha (1216 Å)
    real(wp), parameter :: sigiO_lya  = 0._wp * cs
    real(wp), parameter :: sigiN2_lya = 0._wp * cs
    real(wp), parameter :: sigiO2_lya = 0._wp * cs
    
    !Absorption cross sections (for optical depth tau)
    ! He II
    real(wp), parameter :: sigaO_heii  = sigiO_heii
    real(wp), parameter :: sigaN2_heii = sigiN2_heii
    real(wp), parameter :: sigaO2_heii = sigiO2_heii

    ! He I
    real(wp), parameter :: sigaO_hei  = sigiO_hei
    real(wp), parameter :: sigaN2_hei = sigiN2_hei
    real(wp), parameter :: sigaO2_hei = sigiO2_hei

    ! Ly-beta (important difference)
    real(wp), parameter :: sigaO_lyb  = 0._wp * cs
    real(wp), parameter :: sigaN2_lyb = 0._wp * cs
    real(wp), parameter :: sigaO2_lyb = 1.6_wp * cs

    ! Ly-alpha (model assumption)
    real(wp), parameter :: sigaO_lya  = 0._wp * cs
    real(wp), parameter :: sigaN2_lya = 0._wp * cs
    real(wp), parameter :: sigaO2_lya = 0._wp * cs

    real(wp), dimension(ll) :: brN2i,brN2di,brO2i,brO2di,pepiO,pepiN2i,pepiN2di,pepiO2i,pepiO2di
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol,nN2col,nO2col
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: phototmp
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),ll) :: Iflux
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),llnight) :: Iflux_night
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3),lsp-1) :: photoionization    !don't need a separate rate for electrons
    real(wp) :: alt_km, sza
    real(wp) :: tau_night
    real(wp), dimension(size(nn,1),size(nn,2),size(nn,3)) :: nOcol_vert ,nN2col_vert ,nO2col_vert
    

    !WAVELENGTH BIN BEGINNING AND END (THIS IDEALLY WOULD BE DATA STATEMENTS OR SOME KIND OF STRUCTURE THAT DOESN'T GET REASSIGNED AT EVERY CALL).  Actually all of these array assignments are static...
    lambda1=[0.05, 0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 65.0, &
        79.8, 79.8, 79.8, 91.3, 91.3, 91.3, 97.5, 98.7, 102.7]*1e-9
    lambda2=[0.4, 0.8, 1.8, 3.2, 7.0, 15.5, 22.4, 29.0, 32.0, 54.0, 65.0, 79.8, 79.8, &
         91.3, 91.3, 91.3, 97.5, 97.5, 97.5, 98.7, 102.7, 105.0]*1e-9

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
        
    ! Nighttime resonant/scattered-lines.  Ly-alpha is retained as a placeholder
    ! with zero ionization/absorption cross sections.
    lambda1_night = [29.0_wp, 54.0_wp, 98.7_wp, 105.0_wp] * 1e-9_wp
    lambda2_night = [32.0_wp, 65.0_wp, 102.7_wp, 121.6_wp] * 1e-9_wp

    sigiO_night  = [sigiO_heii,  sigiO_hei,  sigiO_lyb,  sigiO_lya]
    sigiN2_night = [sigiN2_heii, sigiN2_hei, sigiN2_lyb, sigiN2_lya]
    sigiO2_night = [sigiO2_heii, sigiO2_hei, sigiO2_lyb, sigiO2_lya]

    sigaO_night  = [sigaO_heii,  sigaO_hei,  sigaO_lyb,  sigaO_lya]
    sigaN2_night = [sigaN2_heii, sigaN2_hei, sigaN2_lyb, sigaN2_lya]
    sigaO2_night = [sigaO2_heii, sigaO2_hei, sigaO2_lyb, sigaO2_lya]

    ! Branching ratios for night_lines: He II, He I, Ly-beta, Ly-alpha.
    brN2i_night  = [0.908_wp, 1.0_wp,   0.0_wp,   0.0_wp]
    brN2di_night = [0.093_wp, 0.0_wp,   0.0_wp,   0.0_wp]
    brO2i_night  = [0.649_wp, 0.874_wp, 0.613_wp, 0.0_wp]
    brO2di_night = [0.351_wp, 0.108_wp, 0.0_wp,   0.0_wp]
        
    call compute_column_density(nn(:,:,:,1), chi, x, Tninf, gavg, mn(1), nOcol)
    call compute_column_density(nn(:,:,:,2), chi, x, Tninf, gavg, mn(2), nN2col)
    call compute_column_density(nn(:,:,:,3), chi, x, Tninf, gavg, mn(3), nO2col)
    
    Iflux = 0._wp
    do il = 1, ll
      do ix3 = 1, lx3
        do ix2 = 1, lx2
          do ix1 = 1, lx1
          
            !FIXED by adding separate nighttime arrays
            ! FIXME:  
            ! There is a problem here where ll=23 but Iinf only has 22 array entries (solfluxBCs_mod.f90).  
            !   Apparently it doesn't mess up things most
            !   times but this should be fixed as it could have inintended effects depending on mmemory layout.  I would recommend
            !   removing the nighttime only entries (il=23) from the arrays used above for EUVAC and just storing them in other variables.  
            !   This is better because in many cases there will be file-based solar flux inputs that assume 22 bins and we want
            !   the nighttime ionization code to still work with those (e.g. solarfluxBCS.f90 source file in ./boundary_conditions.
            !   So probably the extra 23rd bin data should just be stored in individual variables that get used in the Qnight
            !   calculation.    
!             if (chi(ix1,ix2,ix3) < chi0 + 1._wp*dchi .or. (.not. cfg%flagnightQ)) then    ! don't limit photoionization unless using Qnight
!               Iflux(ix1,ix2,ix3,il) = Iinf(ix1,ix2,ix3,il) * exp( - &
!                    ( sigmaO(il)  * nOcol(ix1,ix2,ix3)  + &
!                      sigmaN2(il) * nN2col(ix1,ix2,ix3) + &
!                      sigmaO2(il) * nO2col(ix1,ix2,ix3) ) )
!             else
!               Iflux(ix1,ix2,ix3,il) = 0._wp
!             end if
            Fchap = Iinf(ix1,ix2,ix3,il) * exp( -( sigmaO(il)  * nOcol(ix1,ix2,ix3)  + &
                   sigmaN2(il) * nN2col(ix1,ix2,ix3) + &
                   sigmaO2(il) * nO2col(ix1,ix2,ix3) ) )

            if (cfg%flagnightQ) then
              w = 0.5_wp * (1._wp - tanh((chi(ix1,ix2,ix3) - chi0)/(2._wp*dchi)))
              Iflux(ix1,ix2,ix3,il) = w * Fchap
            else
              Iflux(ix1,ix2,ix3,il) = Fchap
            end if
          end do
        end do
      end do
    end do
   photoionization = 0._wp

   do il=1,ll
     photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
     photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
     photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
     photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
     photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
   end do
    
    ! Only add nighttime terms if flag is ON
   if (cfg%flagnightQ) then
     call compute_column_density_vertical(nn(:,:,:,1), x, nOcol_vert)
     call compute_column_density_vertical(nn(:,:,:,2), x, nN2col_vert)
     call compute_column_density_vertical(nn(:,:,:,3), x, nO2col_vert)

     Iflux_night = 0._wp
 
     do ix3 = 1, lx3
       do ix2 = 1, lx2
         do ix1 = 1, lx1
           
         sza    = chi(ix1,ix2,ix3)
         alt_km = x%alt(ix1,ix2,ix3) / 1000._wp
         w = 0.5_wp * (1._wp - tanh((sza - chi0)/(2._wp*dchi)))
           
           do il = 1, llnight
             tau_night = sigaO_night(il)*nOcol_vert(ix1,ix2,ix3) + &
                         sigaN2_night(il)*nN2col_vert(ix1,ix2,ix3) + &
                         sigaO2_night(il)*nO2col_vert(ix1,ix2,ix3)

             Fnight = get_nightflux(night_line(il), alt_km, sza) * exp(-tau_night)
             Iflux_night(ix1,ix2,ix3,il) = (1._wp - w) * Fnight
           end do
            
         end do
       end do
     end do
 
     do il = 1, llnight
       photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1) * Iflux_night(:,:,:,il) * sigiO_night(il)
       photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2) * Iflux_night(:,:,:,il) * sigiN2_night(il) * brN2i_night(il)
       photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2) * Iflux_night(:,:,:,il) * sigiN2_night(il) * brN2di_night(il)
       photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3) * Iflux_night(:,:,:,il) * sigiO2_night(il) * brO2i_night(il)
       photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3) * Iflux_night(:,:,:,il) * sigiO2_night(il) * brO2di_night(il)
     end do
   end if

  photoionization(:,:,:,2) = 0._wp
  photoionization(:,:,:,6) = 0._wp

  where (photoionization < 0._wp)
    photoionization = 0._wp
  end where

  do isp=1,lsp-1
    phototmp = photoionization(:,:,:,isp)
    where(x%nullpts)
      phototmp = 0._wp
    end where
    photoionization(:,:,:,isp) = phototmp
  end do

        

!     if (.not. cfg%flagnightQ) then
!     !O COLUMN DENSITY
!     H=kB*Tninf/mn(1)/gavg                         !scalar scale height
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H          !a reduced altitude
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
!     !      Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1._wp-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nOcol=nn(:,:,:,1)*H*Chfn
!         
!         
!     !N2 COLUMN DENSITY
!     H=kB*Tninf/mn(2)/gavg     !all of these temp quantities need to be recomputed for eacb neutral species being ionized
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nN2col=nn(:,:,:,2)*H*Chfn
!         
!         
!     !O2 COLUMN DENSITY
!     H=kB*Tninf/mn(3)/gavg
!     bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H
!     y=sqrt(bigX/2._wp)*abs(cos(chi))
!     Chfn=0
!     where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
!     elsewhere
!       Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!     end where
!     nO2col=nn(:,:,:,3)*H*Chfn
!          
!         
!     !Solar flux at each point on the grid, i.e. the attenuated vacuum flux
!     do il=1,ll
!       do ix3=1,x%lx3
!         do ix2=1,x%lx2
!           do ix1=1,x%lx1 
!           Iflux(ix1,ix2,ix3,il)=Iinf(ix1,ix2,ix3,il)*exp(-(sigmaO(il)*nOcol(ix1,ix2,ix3)+sigmaN2(il)*nN2col(ix1,ix2,ix3)+ &
!                     sigmaO2(il)*nO2col(ix1,ix2,ix3)))
!           end do
!         end do
!       end do
!     end do   
!     
!     !PRIMARY AND SECONDARY IONIZATION RATES
!     photoionization=0
!     
!     !direct O+ production
!     do il=1,ll
!       photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,1)*Iflux(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
!     end do
!         
!     !direct NO+
!     photoionization(:,:,:,2) = 0
!         
!         !direct N2+
!     do il=1,ll
!       photoionization(:,:,:,3)=photoionization(:,:,:,3)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
!     end do
!         
!     !dissociative ionization of N2 leading to N+
!     do il=1,ll
!       photoionization(:,:,:,5)=photoionization(:,:,:,5)+nn(:,:,:,2)*Iflux(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
!     end do
!         
!     !direct O2+
!     do il=1,ll
!       photoionization(:,:,:,4)=photoionization(:,:,:,4)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
!     end do
!         
!     !dissociative ionization of O2 leading to O+
!     do il=1,ll
!       photoionization(:,:,:,1)=photoionization(:,:,:,1)+nn(:,:,:,3)*Iflux(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
!     end do
!         
!     !H+ production
!     photoionization(:,:,:,6) = 0
!     
!     else
!       
!     call compute_column_density(nn(:,:,:,1), chi, x, Tninf, gavg, mn(1), nOcol)
!     call compute_column_density(nn(:,:,:,2), chi, x, Tninf, gavg, mn(2), nN2col)
!     call compute_column_density(nn(:,:,:,3), chi, x, Tninf, gavg, mn(3), nO2col)
!     
!     call compute_column_density_vertical(nn(:,:,:,1), x, nOcol_vert)
!     call compute_column_density_vertical(nn(:,:,:,2), x, nN2col_vert)
!     call compute_column_density_vertical(nn(:,:,:,3), x, nO2col_vert)
!     
!     where (chi < pi/2._wp)
!     nOcol  = nOcol
!     nN2col = nN2col
!     nO2col = nO2col
!     elsewhere
!     nOcol  = nOcol_vert
!     nN2col = nN2col_vert
!     nO2col = nO2col_vert
!     end where
!     
! 
!       ! Flux computation
!       
!     do ix3 = 1, lx3
!      do ix2 = 1, lx2
!       do ix1 = 1, lx1
! 
!       alt_km = x%alt(ix1, ix2, ix3) / 1000._wp
!       sza    = chi(ix1, ix2, ix3)
! 
!         do il = 1, ll
! 
! 
!         ! Daytime direct flux
! 
!         if (sza < chi0) then
!           Fday = Iinf(ix1,ix2,ix3,il) * exp( - &
!                ( sigmaO(il)  * nOcol(ix1,ix2,ix3)  + &
!                  sigmaN2(il) * nN2col(ix1,ix2,ix3) + &
!                  sigmaO2(il) * nO2col(ix1,ix2,ix3) ) )
!         else
!           Fday = 0._wp
!         end if
! 
!         ! Nighttime asymptotic flux
!         
!         Fnight = 0._wp
! 
!         if (il == i_heii) then
!           tau_heii = sigaO_heii  * nOcol_vert(ix1,ix2,ix3)  + &
!                      sigaN2_heii * nN2col_vert(ix1,ix2,ix3) + &
!                      sigaO2_heii * nO2col_vert(ix1,ix2,ix3)
!           Fnight = get_nightflux("heii", alt_km, sza) * exp(-tau_heii)
! 
!         else if (il == i_hei) then
!           tau_hei = sigaO_hei  * nOcol_vert(ix1,ix2,ix3)  + &
!                     sigaN2_hei * nN2col_vert(ix1,ix2,ix3) + &
!                     sigaO2_hei * nO2col_vert(ix1,ix2,ix3)
!           Fnight = get_nightflux("hei", alt_km, sza) * exp(-tau_hei)
! 
!         else if (il == i_lyb) then
!           tau_lyb = sigaO_lyb  * nOcol_vert(ix1,ix2,ix3)  + &
!                     sigaN2_lyb * nN2col_vert(ix1,ix2,ix3) + &
!                     sigaO2_lyb * nO2col_vert(ix1,ix2,ix3)
!           Fnight = get_nightflux("lybeta", alt_km, sza) * exp(-tau_lyb)
! 
!         else if (il == i_lya) then
!           tau_lya = sigaO_lya  * nOcol_vert(ix1,ix2,ix3)  + &
!                     sigaN2_lya * nN2col_vert(ix1,ix2,ix3) + &
!                     sigaO2_lya * nO2col_vert(ix1,ix2,ix3)
!           Fnight = get_nightflux("lyalpha", alt_km, sza) * exp(-tau_lya)
!         end if
! 
! 
!         ! Smooth day-night blending
! 
!         w = 0.5_wp * (1._wp - tanh((sza - chi0)/(2._wp*dchi)))
!         
!         Iflux_day(ix1,ix2,ix3,il)   = w * Fday
!         Iflux_night(ix1,ix2,ix3,il) = (1._wp - w) * Fnight
! 
!         Iflux(ix1,ix2,ix3,il) = Iflux_day(ix1,ix2,ix3,il) + Iflux_night(ix1,ix2,ix3,il)
! 
!         end do
!       end do
!      end do
!     end do
! 
!     ! Ionization rate calculation
!     photoionization = 0._wp
! 
!     do il = 1, ll
! 
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux_day(:,:,:,il)*sigmaO(il)*(1 + pepiO(il))
!         photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux_day(:,:,:,il)*sigmaN2(il)*brN2i(il)*(1 + pepiN2i(il))
!         photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux_day(:,:,:,il)*sigmaN2(il)*brN2di(il)*(1 + pepiN2di(il))
!         photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux_day(:,:,:,il)*sigmaO2(il)*brO2i(il)*(1 + pepiO2i(il))
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux_day(:,:,:,il)*sigmaO2(il)*brO2di(il)*(1 + pepiO2di(il))
!     end do
! 
!     ! Night time
!     do il = 1,ll
!     if (il == i_heii) then    
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux_night(:,:,:,il)*sigiO_heii !*(1 + pepiO(il))
!         photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux_night(:,:,:,il)*sigiN2_heii*brN2i(il) !*(1 + pepiN2i(il))
!         photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux_night(:,:,:,il)*sigiN2_heii*brN2di(il) !*(1 + pepiN2di(il))
!         photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux_night(:,:,:,il)*sigiO2_heii*brO2i(il) !*(1 + pepiO2i(il))
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux_night(:,:,:,il)*sigiO2_heii*brO2di(il) !*(1 + pepiO2di(il))
!     elseif (il == i_hei) then
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,1)*Iflux_night(:,:,:,il)*sigiO_hei !*(1 + pepiO(il))
!         photoionization(:,:,:,3) = photoionization(:,:,:,3) + nn(:,:,:,2)*Iflux_night(:,:,:,il)*sigiN2_hei*brN2i(il) !*(1 + pepiN2i(il))
!         photoionization(:,:,:,5) = photoionization(:,:,:,5) + nn(:,:,:,2)*Iflux_night(:,:,:,il)*sigiN2_hei*brN2di(il) !*(1 + pepiN2di(il))
!         photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux_night(:,:,:,il)*sigiO2_hei*brO2i(il) !*(1 + pepiO2i(il))
!         photoionization(:,:,:,1) = photoionization(:,:,:,1) + nn(:,:,:,3)*Iflux_night(:,:,:,il)*sigiO2_hei*brO2di(il) !*(1 + pepiO2di(il))
!     elseif (il == i_lyb) then
!         ! Ly-beta: only O2 has nonzero cross section
!         photoionization(:,:,:,4) = photoionization(:,:,:,4) + nn(:,:,:,3)*Iflux_night(:,:,:,il)*sigiO2_lyb*brO2i(il) !*(1 + pepiO2i(il))
!     elseif (il == i_lya) then
!         ! Ly-alpha: this contributes zero
!     end if
!     end do
!     end if
! 
!     photoionization(:,:,:,2) = 0._wp
!     photoionization(:,:,:,6) = 0._wp
!     
!     
!     !THERE SHOULD BE SOME CODE HERE TO ZERO OUT THE BELOW-GROUND ALTITUDES.
!     where (photoionization < 0)
!       photoionization = 0
!     end where
!     do isp=1,lsp-1
!       phototmp=photoionization(:,:,:,isp)
!     !  where (x%nullpts>0.9 .and. x%nullpts<1.1)
!       where(x%nullpts)
!         phototmp=0
!       end where
!       photoionization(:,:,:,isp) = phototmp
!     end do
!     
    contains
      subroutine compute_column_density(nn_species, chi, x, Tninf, gavg, mass, n_col)
        real(wp), intent(in) :: nn_species(:,:,:), chi(:,:,:), Tninf, gavg, mass
        class(curvmesh), intent(in) :: x
        real(wp), intent(out) :: n_col(:,:,:)
        real(wp) :: H
        real(wp), dimension(size(nn_species,1), size(nn_species,2), size(nn_species,3)) :: bigX, y, Chfn
      
        H = kB*Tninf/mass/gavg                      !scalar scale height
        bigX=(x%alt(1:lx1,1:lx2,1:lx3)+Re)/H          !a reduced altitude
        y = sqrt(bigX/2._wp)*abs(cos(chi))
        
        where (chi<pi/2._wp)    !where does work with array corresponding elements provided they are conformable (e.g. bigX and y in this case)
            !      Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1._wp-erf(y))    !goodness this creates YUGE errors compared to erfc; left here as a lesson learned
            Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*erfc(y)
        elsewhere
            Chfn=sqrt(pi/2._wp*bigX)*exp(y**2)*(1 + erf(y))
!               Chfn=0._wp
        end where
            n_col = nn_species * H * Chfn      
      end subroutine compute_column_density
      
      !DEBUG
      subroutine compute_column_density_vertical(nn_species, x, n_col)
        real(wp), intent(in) :: nn_species(:,:,:)
        class(curvmesh), intent(in) :: x
        real(wp), intent(out) :: n_col(:,:,:)
        integer :: i1,i2,i3
        real(wp) :: dz

        n_col = 0._wp
  
        do i3 = 1, size(nn_species,3)
          do i2 = 1, size(nn_species,2)
            n_col(size(nn_species,1), i2, i3) = 0._wp   ! column above model top neglected
  
          do i1 = size(nn_species,1)-1, 1, -1
            dz = x%alt(i1+1,i2,i3) - x%alt(i1,i2,i3)   ! meters
            n_col(i1,i2,i3) = n_col(i1+1,i2,i3) + 0.5_wp * dz * &
                            ( nn_species(i1,i2,i3) + nn_species(i1+1,i2,i3) )
          end do
          end do
        end do
      end subroutine compute_column_density_vertical

      function get_nightflux(line_type, alt_km, sza_rad) result(F)
        character(len=*), intent(in) :: line_type
        real(wp), intent(in) :: alt_km, sza_rad
        real(wp) :: F, sza_deg, a, b, c, logF, fscale
        sza_deg = sza_rad * 180.0_wp / pi
        sza_deg = min(max(sza_deg, 90._wp), 180._wp)
!         print *, 'Nighttime SZA =', sza_rad * 180.0_wp / pi
        select case (trim(adjustl(line_type)))
        case ('lyalpha')
        a = 9.7e-05_wp; b = -0.0377_wp; c = 12.7900_wp
        fscale = 1.0_wp
        case ('lybeta')
        a = 1.36e-4_wp; b = -0.0499_wp; c = 10.9299_wp
        fscale = 0.85_wp
    !     print *, 'Case Ly-Beta'
        case ('hei')
        a = -2.4e-05_wp; b = -0.0596_wp; c = 13.5588_wp
        fscale = 0.05_wp
    !     print *, 'Case He I'
        case ('heii')
        a = -1.3e-04_wp; b = 0.0287_wp; c = 6.0416_wp
        fscale = 0.05_wp
    !     print *, 'Case He II'
        case default
        F = 0._wp; return
        end select
    
    
        logF = a * sza_deg**2 + b * sza_deg + c
        F = fscale*10.0_wp**logF*1e4 !m-2s-1
    !      print *, 'Flux for', trim(line_type), 'is', F
    !      print *, 'SZA is:', sza_deg
    !     print *, 'a =', a, ', b =', b, ', c =', c  
      end function get_nightflux
  end function photoionization


  pure function ionrate_fang(W0, PhiWmWm2, alt, nn, Tn, g1, flag_fang, diff_num_flux, kappa, bimax_frac, W0_char)
  real(wp), dimension(:,:), intent(in) :: W0,PhiWmWm2
  real(wp), dimension(:,:,:,:), intent(in) :: nn
  real(wp), dimension(:,:,:), intent(in) :: alt,Tn
  integer, intent(in) :: flag_fang, diff_num_flux
  real(wp), intent(in) :: kappa, bimax_frac, W0_char
  real(wp), dimension(:,:,:), intent(in) :: g1
  real(wp) :: W0keV, PhiW, W0_char_keV
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
          W0_char_keV=W0_char/1e3_wp
    
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
          case (0) ! composite spectrum
            Ptot(:,ix2,ix3) = fang2010_spectrum(PhiW, W0keV, Tn(:,ix2,ix3), massden/1000, meanmass*1000, g1(:,ix2,ix3), &
              diff_num_flux, kappa, bimax_frac, W0_char_keV) * 1e6_wp
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
  
  
  pure function eheating(nn,ionrate,ns)
    !------------------------------------------------------------
    !-------COMPUTE SWARTZ AND NISBET, (1973) ELECTRON HEATING RATES.
    !-------ION ARRAYS (EXCEPT FOR RATES) ARE EXPECTED TO INCLUDE
    !-------GHOST CELLS.
    !------------------------------------------------------------
    
    real(wp), dimension(:,:,:,:), intent(in) :: nn
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
