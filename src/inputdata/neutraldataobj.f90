module neutraldataobj

use phys_consts, only: wp
use inputdataobj, only: inputdata

implicit none (type,external)
public

!> abstract type for neutral data input; this is for 2D or 3D and contains only the data and init
!    other procedures must be defined in concrete class
type, extends(inputdata), abstract :: neutraldata
  ! interpolation site pointer aliases always 3D so defined here in based neutral class
  real(wp), dimension(:,:,:), pointer :: 
  real(wp), dimension(:,:,:), pointer :: dnOinow,dnN2inow,dnO2inow,dvn1inow,dvn2inow,dvn3inow, &
                                           dTninow

  ! bindings for deferred procedures
  contains
end type neutraldata

contains

end module neutraldataobj
