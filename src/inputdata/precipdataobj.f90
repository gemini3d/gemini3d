module precipdataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata
use io, only: ...
implicit none (type, external)

type, extends(inputdata) :: precipdata
  ! coordinate for input precipitation data, and storage
  real(wp), dimension(:), pointer :: mlonp,mlatp
  integer, pointer :: llon,llat
  real(wp), dimension(:,:), pointer :: Qp,E0p

  contains
    init=>init_precip
    set_coordsi=>set_coordsi_precip
    load_data=>load_data_precip
end type precipdata

contains
  !> set pointers to appropriate data arrays (taking into account dimensionality of the problem)
  subroutine init_precip(self,sourcedir,x)
    class(precipdata), intent(inout) :: self
    character, dimension(:), intent(in) :: sourcedir
    class(curvmesh), intent(in) :: x

    ! read the simulation size and grid information from the source directory
    call getsimsize()
    call griddata()
    call self%set_sizes(lc1,lc2,lc3, &
                       0, &
                       0,0,0 &
                       2,0,0 &
                       0, &
                       x )
  end subroutine init_precip


end module precipdataobj
