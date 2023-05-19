module neutraldata3Dobj_fclaw_3Dx

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use reader, only: get_simsize3,get_simsize2,get_grid2,get_precip
use timeutils, only: dateinc,date_filename
use grid, only: gridflag
use geomagnetic, only: ECEFspher2ENU

implicit none (type, external)
private
public :: neutraldata3D_fclaw_3Dx


!> type definition for 3D neutral data that will be provided from a parallel model (i.e. one that runs with GEMINI)
type, extends(neutraldata3D_fclaw) :: neutraldata3D_fclaw_3Dx
  ! for dealing with axisymmetric situations
  real(wp), dimension(:,:,:), allocatable :: proj_ehorzp_e1,proj_ehorzp_e2,proj_ehorzp_e3

  contains
    ! overriding procedures
    procedure :: init_storage
    !procedure :: rotate_winds     ! we don't need to override since base class rotation will work fine for 3D

    ! bindings for deferred procedures
    procedure :: set_coordsi=>set_coordsi_neu3D_fclaw

    ! destructor
    final :: destructor
end type neutraldata3D_fclaw_3Dx


contains
  !> initialize arrays for storing object data once the sizes are set
    subroutine init_storage(self)
    class(neutraldata3D_fclaw_3Dx), intent(inout) :: self
    integer :: lc1,lc2,lc3
    integer :: lc1i,lc2i,lc3i
    integer :: l0D
    integer :: l1Dax1,l1Dax2,l1Dax3
    integer :: l2Dax23,l2Dax12,l2Dax13
    integer :: l3D

    ! check sizes are set
    if (.not. self%flagsizes) error stop 'inpudata:init_storage(); must set sizes before allocations...'

    ! local size variables for convenience
    lc1=self%lc1; lc2=self%lc2; lc3=self%lc3;
    lc1i=self%lc1i; lc2i=self%lc2i; lc3i=self%lc3i;
    l0D=self%l0D
    l1Dax1=self%l1Dax1; l1Dax2=self%l1Dax2; l1Dax3=self%l1Dax3;
    l2Dax23=self%l2Dax23; l2Dax12=self%l2Dax12; l2Dax13=self%l2Dax13;
    l3D=self%l3D

    ! NOTE: type extensions are reponsible for zeroing out any arrays they will use...

!    ! input data coordinate arrays (presume plaid)
    allocate(self%coord1(lc1),self%coord2(lc2),self%coord3(lc3))

    ! interpolation site arrays (note these are flat, i.e. rank 1), if one needed to save space by not allocating unused block
    !   could override this procedure...
    allocate(self%coord1i(lc1i*lc2i*lc3i),self%coord2i(lc1i*lc2i*lc3i),self%coord3i(lc1i*lc2i*lc3i))

!    ! No singleton array objects being allocated by this extension; but this doesn't hurt anything so leave in place
!    ! coordinate sites for singleton axes depend on mangling of data
!    if (self%flagdipmesh) then    ! mangle 2,3 sizes
!      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc3i),self%coord3iax3(lc2i))
!    else
!      allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
!    end if
!    allocate(self%coord2iax23(lc2i*lc3i),self%coord3iax23(lc2i*lc3i))
!    allocate(self%coord1iax13(lc1i*lc3i),self%coord3iax13(lc1i*lc3i))
!    allocate(self%coord1iax12(lc1i*lc2i),self%coord2iax12(lc1i*lc2i))

!    ! allocate object arrays for input data at a reference time.  FIXME: do we even need to store this perm. or can be local to
!    ! load_data?
!    allocate(self%data0D(l0D))
!    allocate(self%data1Dax1(lc1,l1Dax1), self%data1Dax2(lc2,l1Dax2), self%data1Dax3(lc3,l1Dax3))
!    allocate(self%data2Dax23(lc2,lc3,l2Dax23), self%data2Dax12(lc1,lc2,l2Dax12), self%data2Dax13(lc1,lc3,l2Dax13))
!    allocate(self%data3D(lc1,lc2,lc3,l3D))

!    ! allocate object arrays for interpolation sites at reference times
!    allocate(self%data0Di(l0D,2))
!    allocate(self%data1Dax1i(lc1i,l1Dax1,2), self%data1Dax2i(lc2i,l1Dax2,2), self%data1Dax3i(lc3i,l1Dax3,2))
!    allocate(self%data2Dax23i(lc2i,lc3i,l2Dax23,2), self%data2Dax12i(lc1i,lc2i,l2Dax12,2), self%data2Dax13i(lc1i,lc3i,l2Dax13,2))
!    allocate(self%data3Di(lc1i,lc2i,lc3i,l3D,2))

    ! allocate object arrays at interpolation sites for current time.  FIXME: do we even need to store permanently?
    allocate(self%data0Dinow(l0D))
    allocate(self%data1Dax1inow(lc1i,l1Dax1), self%data1Dax2inow(lc2i,l1Dax2), self%data1Dax3inow(lc3i,l1Dax3))
    allocate(self%data2Dax23inow(lc2i,lc3i,l2Dax23), self%data2Dax12inow(lc1i,lc2i,l2Dax12), self%data2Dax13inow(lc1i,lc3i,l2Dax13))
    allocate(self%data3Dinow(lc1i,lc2i,lc3i,l3D))

    !allocate(self%coord1i(lc1i*lc2i*lc3i),self%coord2i(lc1i*lc2i*lc3i),self%coord3i(lc1i*lc2i*lc3i))   
    allocate(self%ximat(lc1i,lc2i,lc3i),self%yimat(lc1i,lc2i,lc3i),self%zimat(lc1i,lc2i,lc3i))
    allocate(self%proj_ezp_e1(lc1i,lc2i,lc3i),self%proj_ezp_e2(lc1i,lc2i,lc3i),self%proj_ezp_e3(lc1i,lc2i,lc3i))
    allocate(self%proj_eyp_e1(lc1i,lc2i,lc3i),self%proj_eyp_e2(lc1i,lc2i,lc3i),self%proj_eyp_e3(lc1i,lc2i,lc3i))
    allocate(self%proj_exp_e1(lc1i,lc2i,lc3i),self%proj_exp_e2(lc1i,lc2i,lc3i),self%proj_exp_e3(lc1i,lc2i,lc3i))

    !FIXME: for when axisymmetric rotations need to be done
    !allocate(self%proj_ehorzp_e1(lc1i,lc2i,lc3i),self%proj_ehorzp_e2(lc1i,lc2i,lc3i),self%proj_ehorzp_e3(lc1i,lc2i,lc3i))

    self%flagalloc=.true.
  end subroutine init_storage



  !> set coordinates for target interpolation points; for neutral inputs we are forced to do some of the property array allocations here
  subroutine set_coordsi_neu3D_fclaw(self,cfg,x)
    class(neutraldata3D_fclaw_3Dx), intent(inout) :: self
    type(gemini_cfg), intent(in) :: cfg
    class(curvmesh), intent(in) :: x
    real(wp) :: theta1,phi1,theta2,phi2,gammarads,theta3,phi3,gamma1,gamma2,phip
    real(wp) :: xp,yp
    real(wp), dimension(3) :: ezp,eyp,tmpvec,exprm
    real(wp), dimension(3) :: erhop    ! FIXME: axisymmetric
    real(wp) :: tmpsca
    integer :: ix1,ix2,ix3,iyn,izn,ixn,iid

    ! Space for coordinate sites and projections in neutraldata3D object
    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i;     ! alias coordinates of interpolation sites

    !Neutral source locations specified in input file, here referenced by spherical magnetic coordinates.
    phi1=cfg%sourcemlon*pi/180
    theta1=pi/2 - cfg%sourcemlat*pi/180

    ! coordinate arrays (ENU)
    call ECEFspher2ENU(x%alt(1:x%lx1,1:x%lx2,1:x%lx3),x%theta(1:x%lx1,1:x%lx2,1:x%lx3),x%phi(1:x%lx1,1:x%lx2,1:x%lx3), &
                          theta1,phi1, &
                          self%ximat,self%yimat,self%zimat)

    do ix3=1,x%lx3
      do ix2=1,x%lx2
        do ix1=1,x%lx1
          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
          ezp=x%er(ix1,ix2,ix3,:)

          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca

          eyp= -x%etheta(ix1,ix2,ix3,:)

          tmpvec=eyp*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=eyp*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_eyp_e3(ix1,ix2,ix3)=tmpsca

          exprm=x%ephi(ix1,ix2,ix3,:)   !for 3D interpolation need to have a unit vector/projection onto x-direction (longitude)

          tmpvec=exprm*x%e1(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e1(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e2(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e2(ix1,ix2,ix3)=tmpsca

          tmpvec=exprm*x%e3(ix1,ix2,ix3,:)
          tmpsca=sum(tmpvec)
          self%proj_exp_e3(ix1,ix2,ix3)=tmpsca


!          !PROJECTIONS FROM NEUTURAL GRID VECTORS TO PLASMA GRID VECTORS
!          !projection factors for mapping from axisymmetric to dipole (go ahead and compute projections so we don't have to do it repeatedly as sim runs
!          ezp=x%er(ix1,ix2,ix3,:)
!          tmpvec=ezp*x%e2(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ezp_e2(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=ezp*x%e1(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ezp_e1(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=ezp*x%e3(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)    !should be zero, but leave it general for now
!          self%proj_ezp_e3(ix1,ix2,ix3)=tmpsca
!
!          erhop=cos(phip)*x%e3(ix1,ix2,ix3,:) - sin(phip)*x%etheta(ix1,ix2,ix3,:)     !unit vector for azimuth (referenced from epicenter - not geocenter!!!) in cartesian geocentric-geomagnetic coords.
!
!          tmpvec=erhop*x%e1(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ehorzp_e1(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=erhop*x%e2(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ehorzp_e2(ix1,ix2,ix3)=tmpsca
!
!          tmpvec=erhop*x%e3(ix1,ix2,ix3,:)
!          tmpsca=sum(tmpvec)
!          self%proj_ehorzp_e3(ix1,ix2,ix3)=tmpsca
        end do
      end do
    end do

    !Assign values for flat lists of grid points
    self%zi=pack(self%zimat,.true.)     !create a flat list of grid points to be used by interpolation functions
    self%yi=pack(self%yimat,.true.)
    self%xi=pack(self%ximat,.true.)

    self%flagcoordsi=.true.
  end subroutine set_coordsi_neu3D_fclaw


!  ! FIXME: overridden with axisymmetric rotation code for now...
!  !> This subroutine takes winds stored in self%dvn?inow and applies a rotational transformation onto the
!  !      grid object for this simulation.  Provided that the horizontal projections have been computed
!  !      correctly the same rotation can be used for axisymmetric and cartesian.
!  subroutine rotate_winds(self)
!    class(neutraldata3D_fclaw_3Dx), intent(inout) :: self
!    integer :: ix1,ix2,ix3
!    real(wp) :: vnhorz,vnz,Tn
!
!    ! do rotations one grid point at a time to cut down on temp storage needed.  Note that until this point there
!    !   shoudl be only zero data stored in vn3 since this class is for 2D data input, instead temperature
!    !   gets stored in the dvn3i variables.
!    do ix3=1,self%lc3i
!      do ix2=1,self%lc2i
!        do ix1=1,self%lc1i
!          vnz=self%dvn1inow(ix1,ix2,ix3)
!          vnhorz=self%dvn2inow(ix1,ix2,ix3)
!          Tn=self%dvn3inow(ix1,ix2,ix3)    ! need to save because it will get overwritten in rotation
!          self%dvn1inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e1(ix1,ix2,ix3) + &
!                                        vnhorz*self%proj_ehorzp_e1(ix1,ix2,ix3)
!          self%dvn2inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e2(ix1,ix2,ix3) + &
!                                        vnhorz*self%proj_ehorzp_e2(ix1,ix2,ix3)
!          self%dvn3inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e3(ix1,ix2,ix3) + &
!                                        vnhorz*self%proj_ehorzp_e3(ix1,ix2,ix3)
!          self%dTninow(ix1,ix2,ix3)=Tn     ! assign saved temperature into correct slot in "output" variables
!        end do
!      end do
!    end do
!  end subroutine rotate_winds


  !> destructor for when object goes out of scope
  subroutine destructor(self)
    type(neutraldata3D_fclaw_3Dx) :: self

    ! deallocate arrays from base inputdata class
    !call self%dissociate_pointers()

    ! null pointers specific to parent neutraldata class
    !call self%dissociate_neutral_pointers()

!    ! I don't know why this causes a segfault...
!    if (associated(self%zlocsi)) deallocate(self%zlocsi)
!    if (associated(self%xlocsi)) deallocate(self%xlocsi)
!    if (associated(self%ylocsi)) deallocate(self%ylocsi)
!    if (associated(self%ilocsi)) deallocate(self%ilocsi)
!    if (associated(self%dataxyzinow)) deallocate(self%dataxyzinow)

    ! due to the nature of this object we cannot rely on base class deallocation
    deallocate(self%data0Dinow)
    deallocate(self%data1Dax1inow, self%data1Dax2inow, self%data1Dax3inow)
    deallocate(self%data2Dax23inow, self%data2Dax12inow, self%data2Dax13inow)
    deallocate(self%data3Dinow)
    deallocate(self%coord1i,self%coord2i,self%coord3i)   

    ! now deallocate arrays specific to this extension
    deallocate(self%proj_ezp_e1,self%proj_ezp_e2,self%proj_ezp_e3)
    deallocate(self%proj_eyp_e1,self%proj_eyp_e2,self%proj_eyp_e3)
    deallocate(self%proj_exp_e1,self%proj_exp_e2,self%proj_exp_e3)
    deallocate(self%ximat,self%yimat,self%zimat)

    ! FIXME: axisymmetric
    !deallocate(self%proj_ehorzp_e1,self%proj_ehorzp_e2,self%proj_ehorzp_e3)   

    ! root has some extra data
!    if (mpi_cfg%myid==0) then
!      deallocate(self%extents,self%indx,self%slabsizes)
!      deallocate(self%xnall,self%ynall)
!    end if

    ! set pointers to null
    nullify(self%xi,self%yi,self%zi);
    !nullify(self%xn,self%yn,self%zn);
    !nullify(self%dnO,self%dnN2,self%dnO2,self%dvnz,self%dvnx,self%dvny,self%dTn)

    self%flagalloc=.false.
    self%flagprimed=.false.
    self%flagcoordsi=.false.
  end subroutine destructor
end module neutraldata3Dobj_fclaw_3Dx
