module neutraldata3Dobj

use phys_consts, only: wp
use neutraldataobj, only: neutraldata
implicit none (type,external)

type, extends(neutraldata) :: neutraldata3D
  ! bindings for deferred procedures
  contains
    init=>init_neu3D
    set_coordsi=>set_coordsi_neu3D
    load_data=>load_data_neu3D
end type neutraldata3D

contains
  !> top set for setting input coords, etc.
  subroutine init_neu3D(self)
    class(neutraldata3D), intent(inout) :: self

    ! sets up data arrays, except for input neutral data coordinates
    call self%initprep_neu()

    ! now neutral data coordinates
    
  end subroutine init_neu3D

end module neutraldata3Dobj
