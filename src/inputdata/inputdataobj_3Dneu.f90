module inputdataobj_3Dneu

use phys_consts, only: wp
use inputdataobj, only: inputdata

implicit none (type,external)
public

type, extends(inputdata) :: inputdata_3Dneu
  ! input coordinate array aliases
  real(wp), dimension(:) :: zn,xn,yn
  ! input data pointer aliases
  real(wp), dimension(:,:,:), pointer :: dnO,dnN2,dnO2,dvny,dvnz,dvnx,dTn
  ! interpolation site pointer aliases
  real(wp), ...

  ! bind deferred procedures
  contains
    procedure :: init=>init_3Dneu
    procedure :: set_coordsi=>set_coordsi_3Dneu
    procedure :: load_data=>load_data_3Dneu
end type inputdata_3Dneu

contains
  subroutine init_3Dneu(self)
    class(inputdata_3Dneu), intent(inout) :: self

    ! alias coordinate arrays
    self%zn=>self%coord1
    self%xn=>self%coord2
    self%yn=>self%coord3

    ! set aliases to point to correct (contiguous) pieces of generic data arrays
    self%dnO=>self%data3D(:,:,:,1)
    self%dnN2=>self%data3D(:,:,:,2)

  end subroutine init_3Dneu(self)

end module inputdataobj_3Dneu
