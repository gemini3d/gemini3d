add_library(advec OBJECT advec_mpi.f90)
target_link_libraries(advec PRIVATE const mpimod grid meshobj)
