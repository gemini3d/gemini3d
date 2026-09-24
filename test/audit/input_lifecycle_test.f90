program input_lifecycle_test
use, intrinsic :: iso_c_binding, only: c_ptr,c_null_ptr,c_associated
use phys_consts, only: wp
use inputdataobj, only: inputdata
use precipdataobj, only: precipdata
use solfluxdataobj, only: solfluxdata
use efielddataobj, only: efielddata
use neutraldataBGobj, only: neutraldataBG
use neutraldataobj, only: neutraldata
use neutraldata2Dobj, only: neutraldata2D
use neutraldata2Dcartobj, only: neutraldata2Dcart
use neutraldata2Daxisymmobj, only: neutraldata2Daxisymm
use neutraldata3Dobj, only: neutraldata3D
use neutraldata3Dobj_mpi, only: neutraldata3D_mpi
use neutraldata3Dobj_geom_mpi, only: neutraldata3D_geom
use neutraldata3Dobj_geog_mpi, only: neutraldata3D_geog
use neutraldata3Dobj_fclaw, only: neutraldata3D_fclaw
use neutraldata3Dobj_fclaw_3Dx, only: neutraldata3D_fclaw_3Dx
use neutraldata3Dobj_fclaw_axisymm, only: neutraldata3D_fclaw_axisymm
use meshobj, only: curvmesh
use meshobj_cart, only: cartmesh
use meshobj_dipole, only: dipolemesh
use gemini3d, only: gemini_grid_dealloc
implicit none
class(inputdata), allocatable :: driver
class(curvmesh), pointer :: mesh
type(c_ptr) :: handle
integer :: iteration,kind,state,xtype
real(wp) :: coords(7)=[1._wp,2._wp,3._wp,4._wp,5._wp,6._wp,7._wp]

do iteration=1,20
  do kind=1,10
    ! Empty, partially constructed, fully allocated, and explicitly released.
    do state=0,3
      select case (kind)
      case (1)
        allocate(precipdata::driver)
      case (2)
        allocate(solfluxdata::driver)
      case (3)
        allocate(efielddata::driver)
      case (4)
        allocate(neutraldataBG::driver)
      case (5)
        allocate(neutraldata2Dcart::driver)
      case (6)
        allocate(neutraldata2Daxisymm::driver)
      case (7)
        allocate(neutraldata3D_geom::driver)
      case (8)
        allocate(neutraldata3D_geog::driver)
      case (9)
        allocate(neutraldata3D_fclaw_3Dx::driver)
      case (10)
        allocate(neutraldata3D_fclaw_axisymm::driver)
      end select
      if (state==1) then
        allocate(driver%lc1,driver%coord1(2),driver%coord1i(8),driver%data3Dinow(2,2,2,7))
      else if (state>=2) then
        call initialize(driver)
      end if
      if (state==3) then
        call driver%dissociate_pointers()
        call assert_released(driver)
        call driver%dissociate_pointers()
      end if
      deallocate(driver)
    end do
  end do

  do kind=1,2
    if (kind==1) then
      allocate(cartmesh::mesh)
    else
      allocate(dipolemesh::mesh)
    end if
    call mesh%set_coords(coords,coords,coords,coords,coords)
    call mesh%init()
    call mesh%init_storage_root()
    if (mod(iteration,2)==0) then
      handle=c_null_ptr; xtype=kind
      call gemini_grid_dealloc(mesh,xtype,handle)
      if (associated(mesh) .or. c_associated(handle)) error stop 'mesh handle not cleared'
    else
      deallocate(mesh)
    end if
  end do
end do
print *, '800 input object and 40 mesh lifecycles passed'

contains

subroutine initialize(self)
  class(inputdata), intent(inout) :: self
  integer :: i

  allocate(self%lc1,self%lc2,self%lc3)
  self%lc1=2; self%lc2=2; self%lc3=2
  self%lc1i=2; self%lc2i=2; self%lc3i=2
  self%l0D=1; self%l1Dax1=1; self%l1Dax2=1; self%l1Dax3=1
  self%l2Dax23=2; self%l2Dax12=1; self%l2Dax13=1; self%l3D=9
  self%flagsizes=.true.; self%flagdatasize=.true.
  select type (self)
  class is (neutraldata3D_fclaw)
    ! Coupled objects allocate only current-time data in their override.
  class is (neutraldata3D)
    allocate(self%coord1(2),self%coord2(2),self%coord3(2))
    allocate(self%coord1i(8),self%coord2i(8),self%coord3i(8))
  class is (neutraldata2D)
    allocate(self%coord1(2),self%coord2(2),self%coord3(2))
    allocate(self%coord1i(8),self%coord2i(8),self%coord3i(8))
  end select
  call self%init_storage()
  do i=1,size(self%coverage)
    allocate(self%coverage(i)%valid(8),source=.true.)
  end do

  select type (self)
  type is (precipdata)
    self%llon=>self%lc2; self%llat=>self%lc3
    self%mlonp=>self%coord2; self%mlatp=>self%coord3
    self%Qp=>self%data2Dax23(:,:,1); self%Qpinow=>self%data2Dax23inow(:,:,1)
  type is (solfluxdata)
    self%llon=>self%lc2; self%llat=>self%lc3; self%lalt=>self%lc1
    self%Iinfinow=>self%data3Dinow
  type is (efielddata)
    self%llon=>self%lc2; self%llat=>self%lc3
    self%flagdirich=>self%data0D(1)
    self%E0xinow=>self%data2Dax23inow(:,:,1)
  type is (neutraldataBG)
    self%llon=>self%lc2; self%llat=>self%lc3; self%lalt=>self%lc1
    self%natmp=>self%data3D; self%natminow=>self%data3Dinow
    allocate(self%proj_ezp_e1(2,2,2))
    self%flagcoordsi=.true. ! remaining allocatables can be absent during teardown
  class is (neutraldata2D)
    self%lzn=>self%lc1; self%lhorzn=>self%lc2; self%lxn=>self%lc3
    self%zn=>self%coord1; self%horzn=>self%coord2
    self%zi=>self%coord1i; self%horzi=>self%coord2i
    allocate(self%proj_ezp_e1(2,2,2),self%zimat(2,2,2))
  class is (neutraldata3D)
    self%lzn=>self%lc1; self%lxn=>self%lc2; self%lyn=>self%lc3
    self%zn=>self%coord1; self%xn=>self%coord2; self%yn=>self%coord3
    self%zi=>self%coord1i; self%xi=>self%coord2i; self%yi=>self%coord3i
  end select
  select type (self)
  class is (neutraldata)
    self%dnOinow=>self%data3Dinow(:,:,:,1)
    self%dvn1inow=>self%data3Dinow(:,:,:,4)
  end select
  select type (self)
  class is (neutraldata3D_mpi)
    allocate(self%extents(2,6),self%indx(2,6),self%slabsizes(2,3))
    allocate(self%xnall(2),self%ynall(2),self%proj_ezp_e1(2,2,2),self%ximat(2,2,2))
  class is (neutraldata3D_fclaw)
    ! Pending borrowed exchange buffers are still owned by the input object.
    allocate(self%zlocsi(8),self%xlocsi(8),self%ylocsi(8),self%ilocsi(8,3),self%dataxyzinow(8,7))
  end select
end subroutine initialize

subroutine assert_released(self)
  class(inputdata), intent(in) :: self
  integer :: i

  if (associated(self%lc1) .or. associated(self%lc2) .or. associated(self%lc3)) error stop 'dimension leak'
  if (associated(self%coord1) .or. associated(self%coord2) .or. associated(self%coord3)) error stop 'source coordinate leak'
  if (associated(self%coord1i) .or. associated(self%coord2i) .or. associated(self%coord3i)) error stop 'target coordinate leak'
  if (associated(self%coord1iax1) .or. associated(self%coord2iax2) .or. associated(self%coord3iax3)) &
    error stop 'axis coordinate leak'
  if (associated(self%coord2iax23) .or. associated(self%coord3iax23)) error stop '23 coordinate leak'
  if (associated(self%coord1iax12) .or. associated(self%coord2iax12)) error stop '12 coordinate leak'
  if (associated(self%coord1iax13) .or. associated(self%coord3iax13)) error stop '13 coordinate leak'
  if (associated(self%data0D) .or. associated(self%data0Di) .or. associated(self%data0Dinow)) error stop 'scalar leak'
  if (associated(self%data1Dax1) .or. associated(self%data1Dax2) .or. associated(self%data1Dax3)) error stop '1D source leak'
  if (associated(self%data1Dax1i) .or. associated(self%data1Dax2i) .or. associated(self%data1Dax3i)) error stop '1D frame leak'
  if (associated(self%data1Dax1inow) .or. associated(self%data1Dax2inow) .or. associated(self%data1Dax3inow)) &
    error stop '1D now leak'
  if (associated(self%data2Dax23) .or. associated(self%data2Dax12) .or. associated(self%data2Dax13)) error stop '2D source leak'
  if (associated(self%data2Dax23i) .or. associated(self%data2Dax12i) .or. associated(self%data2Dax13i)) error stop '2D frame leak'
  if (associated(self%data2Dax23inow) .or. associated(self%data2Dax12inow) .or. associated(self%data2Dax13inow)) &
    error stop '2D now leak'
  if (associated(self%data3D) .or. associated(self%data3Di) .or. associated(self%data3Dinow)) error stop '3D data leak'
  do i=1,size(self%coverage)
    if (allocated(self%coverage(i)%valid)) error stop 'coverage leak'
  end do
  if (self%flagalloc .or. self%flagsizes .or. self%flagdatasize .or. self%flagcoordsi .or. self%flagprimed) &
    error stop 'lifecycle flags not reset'
end subroutine assert_released
end program input_lifecycle_test
