module mumps_rl

implicit none (type, external)
private
public :: @arith@mumps, @arith@mumps_struc

external :: @arith@mumps

include '@arith@mumps_struc.h'

end module mumps_rl

module mumps_interface
use mumps_rl, only : mumps_struc=>@arith@mumps_struc, mumps_exec=>@arith@mumps
implicit none
private
public :: mumps_exec, mumps_struc
end module mumps_interface
