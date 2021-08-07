module neutraldataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata

implicit none (type,external)
public

!> def for 3D neutral data input
type, extends(inputdata), abstract :: neutraldata
  ! input coordinate array aliases, must be defined in extensions
  ! input data pointer aliases
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvny,dvnz,dvnx,dTn
  ! interpolation site pointer aliases
  real(wp), dimension(:,:,:), pointer ::

  ! bindings for deferred procedures
  contains
    procedure :: init=>init_neu
    ! set_coordsi, load_data need to be set by extensions
    !procedure :: set_coordsi=>set_coordsi_3Dneu
    !procedure :: load_data=>load_data_3Dneu
end type neutraldata

contains
  subroutine init_neu(self,sourcedir,dtneu)
    class(neutraldata), intent(inout) :: self
    character, dimension(:), intent(in) :: sourcedir
    integer :: lc1,lc2,lc3
    
    ! read in neutral grid size from file
    call getsimsize(sourcedir,lc1,lc2,lc3)

    ! set sizes
    call self%set_sizes(lc1,lc2,lc3, &
             0, &          ! number scalar parts to dataset
             0, 0, 0, &    ! number 1D data along each axis
             0, 0, 0, &    ! number 2D data
             7)            ! number 3D data

    ! allocate space for arrays
    call self%init_storage()

    ! set the cadence of the dataset input
    self%set_cadence(dtneu)

    ! alias coordinate arrays for internal calculations (as needed)
    !self%zn=>self%coord1
    !self%xn=>self%coord2
    !self%yn=>self%coord3

    ! set aliases to point to correct (contiguous) pieces of generic data arrays
    self%dnO=>self%data3D(:,:,:,1)
    self%dnN2=>self%data3D(:,:,:,2)
    self%dnO2=>self%data3D(:,:,:,3)
    self%dvnz=>self%data3D(:,:,:,4)
    self%dvnx=>self%data3D(:,:,:,5)
    self%dvny=>self%data3D(:,:,:,6)
    self%dTn=>self%data3D(:,:,:,7)
  end subroutine init_neu(self)

end module neutraldataobj
