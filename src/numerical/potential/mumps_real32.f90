module mumps_rl

implicit none (external)
private
public :: smumps, smumps_struc

external :: smumps

include 'smumps_struc.h'

end module mumps_rl

module mumps_interface
use mumps_rl, only : mumps_struc=>smumps_struc, mumps_exec=>smumps
implicit none
private
public :: mumps_exec, mumps_struc
end module mumps_interface
