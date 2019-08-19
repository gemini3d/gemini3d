add_library(potential potential_comm_mumps.f90 potential_mumps.F90 boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PRIVATE calculus collision const grid fsutils interp timeutils
  ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
target_compile_options(potential PRIVATE ${FFLAGS})
target_include_directories(potential PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_definitions(potential PRIVATE REALBITS=${realbits})
set_target_properties(potential PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})