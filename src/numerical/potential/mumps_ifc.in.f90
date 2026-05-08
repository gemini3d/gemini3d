module mumps_rl

implicit none (type, external)
private
public :: @gemini3d_arith@mumps, @gemini3d_arith@mumps_struc

external :: @gemini3d_arith@mumps

include '@gemini3d_arith@mumps_struc.h'

end module mumps_rl

module mumps_interface
use mumps_rl, only : mumps_struc=>@gemini3d_arith@mumps_struc, mumps_exec=>@gemini3d_arith@mumps
implicit none
private
public :: mumps_exec, mumps_struc
end module mumps_interface
