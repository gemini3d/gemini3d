add_library(neutral neutral.f90 atmos.f90)
target_link_libraries(neutral PRIVATE const grid interp mpimod msis timeutils MPI::MPI_Fortran)
target_include_directories(neutral INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(neutral PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)