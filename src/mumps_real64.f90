module mumps_rl
implicit none
private
public :: dmumps, dmumps_struc

external :: dmumps

include 'dmumps_struc.h'

end module mumps_rl

module mumps_interface
use mumps_rl, only : mumps_struc=>dmumps_struc, mumps_exec=>dmumps
implicit none
private
public :: mumps_exec, mumps_struc
end module mumps_interface
