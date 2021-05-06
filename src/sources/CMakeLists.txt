add_library(sources OBJECT sources.f90 sources_mpi.f90)
target_link_libraries(sources PRIVATE calculus collision const grid meshobj mpimod)
