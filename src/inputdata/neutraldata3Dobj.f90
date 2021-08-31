module neutraldata3Dobj

use phys_consts, only: wp
use neutraldataobj, only: neutraldata
implicit none (type,external)

!> type definition for 3D neutral data
type, extends(neutraldata) :: neutraldata3D
  ! source data coordinate pointers
  real(wp), dimension(:), pointer :: xn,yn,zn

  ! source data pointers
  real(wp), dimension(:,:,:), pointer :: 
  contains
    ! bindings for deferred procedures
    procedure :: init=>init_neu3D
    procedure :: set_coordsi=>set_coordsi_neu3D
    procedure :: load_data=>load_data_neu3D
    procedure :: load_grid=>load_grid_neu3D
    procedure :: load_size=>load_size_neu3D 
end type neutraldata3D

contains
  !> initialize storage for this type of neutral input data
  subroutine init_neu3D(self,sourcedir,dtneu)
    class(neutraldata3D), intent(inout) :: self
    character, dimension(:), intent(in) :: sourcedir
    integer :: lc1,lc2,lc3
    
    ! read in neutral grid size from file (io module)
    call getsimsize(sourcedir,lc1,lc2,lc3)

    ! set sizes, we have 7 arrays all 3D (irrespective of 2D vs. 3D neutral input)
    call self%set_sizes(lc1,lc2,lc3, &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data
             7)            ! number 3D data

    ! allocate space for arrays
    call self%init_storage()
    call self%set_cadence(dtneu)

    ! alias coordinate arrays for internal calculations (as needed)
    self%zn=>self%coord1; self%xn=>self%coord2; self%yn=>self%coord3;

    ! set aliases to point to correct source data arrays
    self%dnO=>self%data3D(:,:,:,1)
    self%dnN2=>self%data3D(:,:,:,2)
    self%dnO2=>self%data3D(:,:,:,3)
    self%dvnz=>self%data3D(:,:,:,4)
    self%dvnx=>self%data3D(:,:,:,5)
    self%dvny=>self%data3D(:,:,:,6)
    self%dTn=>self%data3D(:,:,:,7)


    ! set aliases for interpolated data that is "outward facing"
    self%dnOinow=>self%data3Dinow(:,:,:,1)
    self%dnN2inow=>self%data3Dinow(:,:,:,2)
    self%dnO2inow=>self%data3Dinow(:,:,:,3)
    self%dvn1inow=>self%data3Dinow(:,:,:,4)
    self%dvn2inow=>self%data3Dinow(:,:,:,5)
    self%dvn3inow=>self%data3Dinow(:,:,:,6)
    self%dTninow=>self%data3Dinow(:,:,:,7)
  end subroutine init_neu3D
end module neutraldata3Dobj
