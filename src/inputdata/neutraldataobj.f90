module neutraldataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata

implicit none (type,external)
private
public :: neutraldata

!> abstract type for neutral data input; this is for 2D or 3D and contains only the data and init
!    other procedures must be defined in concrete class
type, extends(inputdata), abstract :: neutraldata
  ! interpolation site pointer aliases always 3D so defined here in based neutral class
  real(wp), dimension(:,:,:), pointer ::  dnOiprev,dnN2iprev,dnO2iprev,dvn1iprev,dvn2iprev,dvn3iprev, &
                                           dTniprev
  real(wp), dimension(:,:,:), pointer ::  dnOinext,dnN2inext,dnO2inext,dvn1inext,dvn2inext,dvn3inext, &
                                           dTninext
  real(wp), dimension(:,:,:), pointer :: dnOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow, &
                                           dTninow
  contains
    procedure :: setptrs_grid
    procedure :: dissociate_neutral_pointers
end type neutraldata

contains
  !> set pointer variables to locations for storage of interpolated data (3D always for neutral input)
  subroutine setptrs_grid(self)
    class(neutraldata), intent(inout) :: self

    ! set aliases for prev data
    self%dnOiprev=>self%data3Di(:,:,:,1,1)
    self%dnN2iprev=>self%data3Di(:,:,:,2,1)
    self%dnO2iprev=>self%data3Di(:,:,:,3,1)
    self%dvn1iprev=>self%data3Di(:,:,:,4,1)
    self%dvn2iprev=>self%data3Di(:,:,:,5,1)
    self%dvn3iprev=>self%data3Di(:,:,:,6,1)
    self%dTniprev=>self%data3Di(:,:,:,7,1)

    ! set pointers for next data
    self%dnOinext=>self%data3Di(:,:,:,1,2)
    self%dnN2inext=>self%data3Di(:,:,:,2,2)
    self%dnO2inext=>self%data3Di(:,:,:,3,2)
    self%dvn1inext=>self%data3Di(:,:,:,4,2)
    self%dvn2inext=>self%data3Di(:,:,:,5,2)
    self%dvn3inext=>self%data3Di(:,:,:,6,2)
    self%dTninext=>self%data3Di(:,:,:,7,2)

    ! set aliases for interpolated data that is "outward facing"
    self%dnOinow=>self%data3Dinow(:,:,:,1)
    self%dnN2inow=>self%data3Dinow(:,:,:,2)
    self%dnO2inow=>self%data3Dinow(:,:,:,3)
    self%dvn1inow=>self%data3Dinow(:,:,:,4)
    self%dvn2inow=>self%data3Dinow(:,:,:,5)
    self%dvn3inow=>self%data3Dinow(:,:,:,6)
    self%dTninow=>self%data3Dinow(:,:,:,7)
  end subroutine setptrs_grid


  !> nullify neutral pointers (dealloc should occur from base class); unclear whether fortran standard automatically
  !    calls for setting pointers to null vs. undefined.
  subroutine dissociate_neutral_pointers(self)
    class(neutraldata), intent(inout) :: self

    nullify(self%dnOiprev,self%dnN2iprev,self%dnO2iprev,self%dvn1iprev,self%dvn2iprev,self%dvn3iprev, &
                                           self%dTniprev)
    nullify(self%dnOinext,self%dnN2inext,self%dnO2inext,self%dvn1inext,self%dvn2inext,self%dvn3inext, &
                                           self%dTninext)
    nullify(self%dnOinow,self%dnN2inow,self%dnO2inow,self%dvn1inow,self%dvn2inow,self%dvn3inow, &
                                           self%dTninow)
  end subroutine dissociate_neutral_pointers
end module neutraldataobj
