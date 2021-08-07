module neutraldataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata

implicit none (type,external)
public

!> abstract type for neutral data input; this is for 2D or 3D and contains only the data and init
!    other procedures must be defined in concrete class
type, extends(inputdata), abstract :: neutraldata
  ! input coordinate array aliases, must be defined in extensions
  ! input data pointer aliases
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvnz,dvnx,dvny,dTn
  ! interpolation site pointer aliases
  real(wp), dimension(:,:,:), pointer :: dnOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow, &
                                           dTninow

  ! bindings for deferred procedures
  contains
    procedure :: initprep_neu
    ! deferred procedures need to be set by extensions
    !procedure :: init=>...
    !procedure :: set_coordsi=>set_coordsi_3Dneu
    !procedure :: load_data=>load_data_3Dneu
end type neutraldata

contains
  subroutine initprep_neu(self,sourcedir,dtneu)
    class(neutraldata), intent(inout) :: self
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

    ! set aliases for interpolated data that is "outward facing"
    self%dnOinow=>self%data3Dinow(:,:,:,1)
    self%dnN2inow=>self%data3Dinow(:,:,:,2)
    self%dnO2inow=>self%data3Dinow(:,:,:,3)
    self%dvn1inow=>self%data3Dinow(:,:,:,4)
    self%dvn2inow=>self%data3Dinow(:,:,:,5)
    self%dvn3inow=>self%data3Dinow(:,:,:,6)
    self%dTninow=>self%data3Dinow(:,:,:,7)
  end subroutine initprep_neu(self)
end module neutraldataobj
