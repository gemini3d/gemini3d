add_library(neutral OBJECT neutral.f90 atmos.f90)
target_link_libraries(neutral PRIVATE const reader grid mesh interp mpimod msis timeutils MPI::MPI_Fortran)

# --- for setting up an equilibrium simulation --

add_executable(msis_setup msis_driver.f90)
target_link_libraries(msis_setup PRIVATE msis)
set_target_properties(msis_setup PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})