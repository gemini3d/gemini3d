module solfluxBCs_mod

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

use reader, only: get_simsize2, get_precip, get_grid2
use phys_consts, only: pi,wp, debug
use grid, only : lx1,lx2,lx3
use meshobj, only: curvmesh
use interpolation, only : interp1,interp2
use timeutils, only : dateinc, date_filename, find_lastdate
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi
use gemini3d_config, only: gemini_cfg
use solfluxdataobj, only: solfluxdata

implicit none (type, external)
private
public :: solfluxBCs_fileinput, solfluxBCs, init_solfluxinput

contains
  !> initialize variables to hold input file precipitation information, must be called by all workers at the same time
  subroutine init_solfluxinput(dt,cfg,ymd,UTsec,x,Iinf,solflux)
    real(wp), intent(in) :: dt
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    class(curvmesh), intent(in) :: x
    real(wp), dimension(:,:,:,:), intent(inout) :: Iinf
    type(solfluxdata), intent(inout) :: solflux

    if (cfg%flagsolfluxfile==1) then    !all workers must have this info
      call solflux%init(cfg,cfg%precdir,x,dt,cfg%dtprec,ymd,UTsec)
    end if
  end subroutine init_solfluxinput


  !> get latest file input precipitation information
  subroutine solfluxBCs_fileinput(dtmodel,t,cfg,ymd,UTsec,x,Iinf,solflux)
    real(wp), intent(in) :: dtmodel
    real(wp), intent(in) :: t
    type(gemini_cfg), intent(in) :: cfg
    integer, dimension(3), intent(in) :: ymd
    !! date for which we wish to calculate perturbations
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), intent(inout) :: Iinf
    !! intent(out)
    !! last dimension is the number of particle populations
    class(curvmesh), intent(in) :: x
    type(solfluxdata), intent(inout) :: solflux
    integer :: ix2,ix3

    ! disturbance precipitation from file input
    call solflux%update(cfg,dtmodel,t,x,ymd,UTsec)

    ! set output arrays; note that this is making a copy of the data stored the precipdata object.  We
    !   will assume for now that this doesn't incur too much memory overhead
    Iinf(:,:,:,:)=solflux%Iinfinow(:,:,:,:)
  end subroutine solfluxBCs_fileinput


  !> This is the default subroutine that is called for *vacuum* (large altitude limit) solar flux values above each position on the grid 
  subroutine solfluxBCs(cfg,x,ymd,UTsec,Iinf)
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    integer, dimension(3), intent(in) :: ymd
    real(wp), intent(in) :: UTsec
    real(wp), dimension(:,:,:,:), intent(inout) :: Iinf
    !! intent(out)
    integer, parameter :: ll=22     !number of wavelength bins   
    real(wp), dimension(ll) :: fref,Aeuv   
    !> 2D mask below
    !real(wp), parameter :: ecglat =33.0, ecglong=255.0, ecwidth=5.0, ectime=63000.0, ecdtime=1800.0, maskmax=0.9
    real(wp), dimension(2), parameter :: ecglat=[49.0,-5.0], ecglon=[213,331], ectime=[56700.0,69183.33]
    real(wp), parameter :: ecwidth=5.0, ecdtime=1800.0, maskmax=0.9
    real(wp) :: ecglatnow,ecglonnow
    logical, parameter :: flagmask=.false.    ! hardcoded toggle for eclipse
    integer :: il,ix1,ix2,ix3
    real(wp) :: maskval
    real(wp), dimension(ll) :: Iinfref
    real(wp) :: f107,f107a

    f107=cfg%activ(2)
    f107a=cfg%activ(1)

    !EUVAC FLUX VALUES
    fref=[5.01e1, 1e4, 2e6, 2.85e7, 5.326e8, 1.27e9, 5.612e9, 4.342e9, 8.380e9, &
        2.861e9, 4.83e9, 1.459e9, 1.142e9, 2.364e9, 3.655e9, 8.448e8, 3.818e8, &
        1.028e9, 7.156e8, 4.482e9, 4.419e9, 4.235e9]*1e4        !convert to m^-2 s^-1
    Aeuv=[6.24e-1, 3.71e-1, 2e-1, 6.247e-2, 1.343e-2, 9.182e-3, 1.433e-2, 2.575e-2, &
        7.059e-3, 1.458e-2, 5.857e-3, 5.719e-3, 3.680e-3, 5.310e-3, 5.261e-3, 5.437e-3, &
        4.915e-3, 4.995e-3, 4.422e-3, 3.950e-3, 5.021e-3, 4.825e-3]

    !IRRADIANCE ACCORDING TO [RICHARDS, 1994]
    Iinfref=fref*(1 + Aeuv*(0.5_wp*(f107+f107a)-80._wp))

    !PHOTON FLUX
    do il=1,ll
      ! if applying an eclipse mask we need to loop over the positions and determine the mask value for each before
      !   computing photon fluxes. 
      do ix3=1,x%lx3
        do ix2=1,x%lx2
          do ix1=1,x%lx1 
            if (flagmask) then
              ecglatnow=ecglat(1)+(ecglat(2)-ecglat(1))/(ectime(2)-ectime(1))
              ecglonnow=ecglon(1)+(ecglon(2)-ecglon(1))/(ectime(2)-ectime(1))
    !> 2D mask
    !          maskval=maskmax*exp(-(x%glat(ix1,ix2,ix3)-ecglat)**2/2/ecwidth**2)* &
    !                    !exp(-(x%glon(ix1,ix2,ix3)-ecglong)**2/2/ecwidth**2)* &
    !                    exp(-(UTsec-ectime)**2/2/ecdtime**2)
    
              if (UTsec>ectime(1) .and. UTsec<ectime(2)) then
                maskval=maskmax*exp(-(x%glat(ix1,ix2,ix3)-ecglatnow)**2/2/ecwidth**2)* &
                          exp(-(x%glon(ix1,ix2,ix3)-ecglonnow)**2/2/ecwidth**2)
              else
                maskval=0.0
              end if
            else
              maskval=0.0
            end if
            Iinf(ix1,ix2,ix3,il)=(1.0-maskval)*Iinfref(il)
          end do
        end do
      end do
    end do   
  end subroutine solfluxBCs
end module solfluxBCs_mod
