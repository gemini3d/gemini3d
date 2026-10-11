module diffusion_features
implicit none(type,external)
private
public :: midpoint_conduction_enabled
logical, parameter :: midpoint_conduction_enabled=@midpoint_conduction_fortran@
end module diffusion_features
