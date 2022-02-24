add_library(advec OBJECT advec.f90)
target_link_libraries(advec PRIVATE const grid meshobj)

add_library(advec_mpi OBJECT advec_mpi.f90)
target_link_libraries(advec_mpi PRIVATE const mpimod grid meshobj)
