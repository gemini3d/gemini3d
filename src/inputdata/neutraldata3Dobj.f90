!> This module contains the parent object for all 3D neutral perturbation inputdata objects.
module neutraldata3Dobj

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use phys_consts, only: wp,debug,pi,Re
use inputdataobj, only: inputdata
use neutraldataobj, only: neutraldata
use meshobj, only: curvmesh
use gemini3d_config, only: gemini_cfg
use timeutils, only: dateinc,date_filename

implicit none (type, external)
private
public :: neutraldata3D

type, abstract, extends(neutraldata) :: neutraldata3D
  ! source data coordinate pointers
  real(wp), dimension(:), pointer :: xn,yn,zn
  integer, pointer :: lxn,lyn,lzn
  real(wp), dimension(:), allocatable :: xnall,ynall
  integer :: lxnall,lynall

  ! work arrays needed by various procedures re: target coordinates
  real(wp), dimension(:,:,:), allocatable :: ximat,yimat,zimat
  real(wp), dimension(:), pointer :: zi,xi,yi

  ! source data pointers
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvnz,dvnx,dvny,dTn

  ! projection factors needed to rotate input data onto grid
  real(wp), dimension(:,:,:), allocatable :: proj_ezp_e1,proj_ezp_e2,proj_ezp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_eyp_e1,proj_eyp_e2,proj_eyp_e3
  real(wp), dimension(:,:,:), allocatable :: proj_exp_e1,proj_exp_e2,proj_exp_e3

  contains
    ! unique to this class
    procedure :: rotate_winds
    
    ! overriding procedures
    procedure :: init_storage

    ! deferred procedures (can be overridden by child class as needed)
end type neutraldata3D

contains
  !> create storage for arrays needed specifically for 3D neutral input calculations, overrides the base class procedure
  subroutine init_storage(self)
    class(neutraldata3D), intent(inout) :: self
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

    ! input data coordinate arrays are set by load_gridandsize()

    ! allocate target coords, for neutral3D the standard set (coord1i, etc.) are done in set_coordsi()
    allocate(self%coord1iax1(lc1i),self%coord2iax2(lc2i),self%coord3iax3(lc3i))
    allocate(self%coord2iax23(lc2i*lc3i),self%coord3iax23(lc2i*lc3i))
    allocate(self%coord1iax13(lc1i*lc3i),self%coord3iax13(lc1i*lc3i))
    allocate(self%coord1iax12(lc1i*lc2i),self%coord2iax12(lc1i*lc2i))

    ! allocate object arrays for input data at a reference time.  FIXME: do we even need to store this perm. or can be local to
    ! load_data?
    allocate(self%data0D(l0D))
    allocate(self%data1Dax1(lc1,l1Dax1), self%data1Dax2(lc2,l1Dax2), self%data1Dax3(lc3,l1Dax3))
    allocate(self%data2Dax23(lc2,lc3,l2Dax23), self%data2Dax12(lc1,lc2,l2Dax12), self%data2Dax13(lc1,lc3,l2Dax13))
    allocate(self%data3D(lc1,lc2,lc3,l3D))

    ! allocate object arrays for interpolation sites at reference times
    allocate(self%data0Di(l0D,2))
    allocate(self%data1Dax1i(lc1i,l1Dax1,2), self%data1Dax2i(lc2i,l1Dax2,2), self%data1Dax3i(lc3i,l1Dax3,2))
    allocate(self%data2Dax23i(lc2i,lc3i,l2Dax23,2), self%data2Dax12i(lc1i,lc2i,l2Dax12,2), self%data2Dax13i(lc1i,lc3i,l2Dax13,2))
    allocate(self%data3Di(lc1i,lc2i,lc3i,l3D,2))

    ! allocate object arrays at interpolation sites for current time.  FIXME: do we even need to store permanently?
    allocate(self%data0Dinow(l0D))
    allocate(self%data1Dax1inow(lc1i,l1Dax1), self%data1Dax2inow(lc2i,l1Dax2), self%data1Dax3inow(lc3i,l1Dax3))
    allocate(self%data2Dax23inow(lc2i,lc3i,l2Dax23), self%data2Dax12inow(lc1i,lc2i,l2Dax12), self%data2Dax13inow(lc1i,lc3i,l2Dax13))
    allocate(self%data3Dinow(lc1i,lc2i,lc3i,l3D))

    self%flagalloc=.true.
  end subroutine init_storage


  !> This subroutine takes winds stored in self%dvn?inow and applies a rotational transformation onto the
  !      grid object for this simulation
  subroutine rotate_winds(self)
    class(neutraldata3D), intent(inout) :: self
    integer :: ix1,ix2,ix3
    real(wp) :: vnx,vny,vnz

    ! do rotations one grid point at a time to cut down on temp storage needed
    do ix3=1,self%lc3i
      do ix2=1,self%lc2i
        do ix1=1,self%lc1i
          vnz=self%dvn1inow(ix1,ix2,ix3)
          vnx=self%dvn2inow(ix1,ix2,ix3)
          vny=self%dvn3inow(ix1,ix2,ix3)
          self%dvn1inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e1(ix1,ix2,ix3) + vnx*self%proj_exp_e1(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e1(ix1,ix2,ix3)
          self%dvn2inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e2(ix1,ix2,ix3) + vnx*self%proj_exp_e2(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e2(ix1,ix2,ix3)
          self%dvn3inow(ix1,ix2,ix3)=vnz*self%proj_ezp_e3(ix1,ix2,ix3) + vnx*self%proj_exp_e3(ix1,ix2,ix3) + &
                                        vny*self%proj_eyp_e3(ix1,ix2,ix3)
        end do
      end do
    end do
  end subroutine rotate_winds
end module neutraldata3Dobj
