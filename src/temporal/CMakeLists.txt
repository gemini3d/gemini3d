set_directory_properties(PROPERTIES LABELS temporal)

add_library(timeutils OBJECT timeutils.f90)
target_link_libraries(timeutils PRIVATE const)

add_library(temporal OBJECT temporal.f90)
target_link_libraries(temporal PRIVATE config const meshobj mpimod MPI::MPI_Fortran)
