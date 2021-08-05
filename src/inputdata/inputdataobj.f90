module inputdataobj

use phys_consts, only : wp

implicit none (type, external)
public

!> this is a generic class for an data object being input into the model and interpolated in space and time
type, abstract :: inputdata
  real(wp), dimension(:), pointer :: coord1,coord2,coord3     ! coordinates for the source data (interpolant coords)
  integer :: lc1,lc2,lc3                                      ! dataset length along the 3 coordinate axes

  !! we only store data that have already been spatially interpolated
  real(wp), dimension(:,:), pointer :: data0Di                       ! array for storing a "stack" of scalar data (only interpolated in time)
                                                                       !  last axis is for prev,next copies of data for interp in time
                                                                       !  second to last axis is for number of datasets of this dimension
  integer :: axes1D
  integer :: l0D                                                     ! length/number of scalar datasets
  real(wp), dimension(:,:,:), pointer :: data1Di                     ! array for storing series of 1D data, array varies along non-singleton axis
  integer :: l1D
  real(wp), dimension(:,:,:,:), pointer :: data2Di                   ! array for storing series of 2D data, varies along two non-singleton axes
  integer :: l2D
  real(wp), dimension(:,:,:,:,:), pointer :: data3Di                 ! array for storing series of 3D data
  integer :: l3D
  real(wp), dimension(:,:,:), pointer :: coord1i,coord2i,coord3i     ! coordinates of the interpolation sites
  integer :: lc1i,lc2i,lc3i                                          ! dataset length along the 3 coordinate axes

  real(wp), dimension(2) :: tref                                     ! times for two input frames bracketting current time

  contains
    procedure :: setsizes     ! initiate sizes for coordinate axes and number of datasets of different dimensionality
    procedure :: setcoords    ! fill interpolant coordinate arrays
    procedure :: init_storage ! wrapper routine to set up arrays and sizes
    procedure :: update       ! check to see if new file needs to be read and read accordingly (will need to call deferred loaddata)
    procedure :: timeinterp   ! interpolate in time based on data presently loaded into spatial arrays
    procedure :: spaceinterp  ! interpolate spatially
    procedure :: dissociate_pointers   ! clear out memory and reset and allocation status flags

    procedure(initproc), deferred :: init        ! create storage (sub), and prime data as needed
    procedure(loadproc), deferred :: loaddata    ! read data from file (possibly one array at a time) and spatially interpolate and store it in the appropriate arrays

end type inputdata


!> interfaces for deferred procedures
abstract interface
  subroutine initproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine init
  subroutine loadproc(self)
    import inputdata
    class(inputdata), intent(inout) :: self
  end subroutine loadproc
end interface

contains

end module inputdataobj
