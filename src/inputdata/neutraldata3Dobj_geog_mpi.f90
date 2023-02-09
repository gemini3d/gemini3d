module neutraldata3Dobj_geog_mpi

use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj_mpi, only: neutraldata3D_mpi
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use geomagnetic, only: geomag2geog, ECEFspher2ENU
use mpimod, only: mpi_realprec, mpi_cfg, tag=>gemini_mpi

use mpi_f08, only: mpi_send,mpi_recv,mpi_integer,mpi_comm_world

implicit none (type, external)
private
public :: neutraldata3D_geog

!> type definition for 3D neutral data in geographic coordinates
type, extends(neutraldata3D_mpi) :: neutraldata3D_geog
  !! all data use parent class pointers/arrays
  contains
    !! new deferred binding
    procedure :: set_coordsi=>set_coordsi_neu3D_geog

    !! destructor
    final :: destructor
end type neutraldata3D_geog

contains
  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
  subroutine set_coordsi_neu3D_geog(self,cfg,x)
    class(neutraldata3D_geog), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: theta1,phi1,theta2,phi2,theta3,phi3,gamma1,gamma2
    real(wp) :: xp,yp
    real(wp), dimension(3) :: ezp,eyp,tmpvec,exprm
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3
    real(wp) :: glon1,glat1
    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3,3) :: ealt,eglat,eglon
    real(wp), dimension(1:x%lx1,1:x%lx2,1:x%lx3) :: thetageog,phigeog


    ! Space for coordinate sites and projections in neutraldata3D object
    allocate(self%coord1i(x%lx1*x%lx2*x%lx3),self%coord2i(x%lx1*x%lx2*x%lx3),self%coord3i(x%lx1*x%lx2*x%lx3))
    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i;     ! coordinates of interpolation sites
    allocate(self%ximat(x%lx1,x%lx2,x%lx3),self%yimat(x%lx1,x%lx2,x%lx3),self%zimat(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_ezp_e1(x%lx1,x%lx2,x%lx3),self%proj_ezp_e2(x%lx1,x%lx2,x%lx3),self%proj_ezp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_eyp_e1(x%lx1,x%lx2,x%lx3),self%proj_eyp_e2(x%lx1,x%lx2,x%lx3),self%proj_eyp_e3(x%lx1,x%lx2,x%lx3))
    allocate(self%proj_exp_e1(x%lx1,x%lx2,x%lx3),self%proj_exp_e2(x%lx1,x%lx2,x%lx3),self%proj_exp_e3(x%lx1,x%lx2,x%lx3))

    !Neutral source locations specified in input file, here we must convert geomag. (in input file) to geographic geographic coordinates
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2 - cfg%sourcemlat*pi/180
    call geomag2geog(phi1,theta1,glon1,glat1)
    phi1=glon1*pi/180
    theta1=pi/2 - glat1*pi/180

    !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
    if (mpi_cfg%myid==0) print*, 'Getting unit vectors for geographic directions on mesh...'
    call x%calc_unitvec_geo(ealt,eglon,eglat)

    !Convert plasma simulation grid locations to z,rho values to be used in interoplation.  altitude ~ zi; lat/lon --> rhoi.  Also compute unit vectors and projections
    if (mpi_cfg%myid==0) then
      print *, 'Computing alt,radial distance values for plasma grid and completing rotations, using geographic coordinates...'
    end if
    thetageog=pi/2._wp - x%glat(1:x%lx1,1:x%lx2,1:x%lx3)*pi/180._wp
    phigeog=x%glon(1:x%lx1,1:x%lx2,1:x%lx3)*pi/180._wp
    call ECEFspher2ENU(x%alt(1:x%lx1,1:x%lx2,1:x%lx3),thetageog,phigeog,theta1,phi1,self%ximat,self%yimat,self%zimat)

    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections as well)
          ezp=ealt(ix1,ix2,ix3,:)
          !ezp=x%er(ix1,ix2,ix3,:)

          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca

          ! we now need geographic unit vectors which we can get from our grid methods
          eyp=eglat(ix1,ix2,ix3,:)
          !eyp= -x%etheta(ix1,ix2,ix3,:)

          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca

          exprm=eglon(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)
          !exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)

          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e3(ix1,ix2,ix3)=tmpsca
        end do
      end do
    end do

    !Assign values for flat lists of grid points
    if (mpi_cfg%myid==0) then
      print*, '...Packing interpolation target points...'
    end if
    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
    self%yi=pack(self%yimat,.true.)
    self%xi=pack(self%ximat,.true.)

    ! FIXME: do we need to have the new grid code clear its unit vectors?  Or maybe this isn't a huge waste of memory???
    if (mpi_cfg%myid==0) then
      print*, '...Clearing out unit vectors (after projections)...'
    end if
    !call clear_unitvecs(x)

    if(mpi_cfg%myid==0) then
      print*, 'Interpolation coords:  ',minval(self%zi),maxval(self%zi), &
                                        minval(self%xi),maxval(self%xi), &
                                        minval(self%yi),maxval(self%yi)
      print*, 'Projection checking:  ',minval(self%proj_exp_e1),maxval(self%proj_exp_e1), &
                                       minval(self%proj_exp_e2),maxval(self%proj_exp_e2), &
                                       minval(self%proj_exp_e3),maxval(self%proj_exp_e3)
    end if

    self%flagcoordsi=.true.
  end subroutine set_coordsi_neu3D_geog


  !> destructor for when object goes out of scope
  subroutine destructor(self)
    type(neutraldata3D_geog) :: self

    ! deallocate arrays from base inputdata class
    call self%dissociate_pointers()

    ! null pointers specific to parent neutraldata class
    call self%dissociate_neutral_pointers()

    ! now deallocate arrays specific to this extension
    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
    deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
    deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
    deallocate(self%ximat,self%yimat,self%zimat)

    ! root has some extra data
    if (mpi_cfg%myid==0) then
      deallocate(self%extents,self%indx,self%slabsizes)
      deallocate(self%xnall,self%ynall)
    end if

    ! set pointers to null
    nullify(self%xi,self%yi,self%zi);
    nullify(self%xn,self%yn,self%zn);
    nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnx,self%dvny,self%dTn)
  end subroutine destructor

end module neutraldata3Dobj_geog_mpi
