add_library(sources sources.f90 sources_mpi.f90)
target_link_libraries(sources PRIVATE calculus collision const grid mesh mpimod)
target_include_directories(sources INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(sources PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
