add_library(mpi_c_stub mpic.c)
target_include_directories(mpi_c_stub INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_library(MPI::MPI_C ALIAS mpi_c_stub)

add_library(elapse OBJECT elapse.c)
target_include_directories(elapse INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(mpi_f_stub mpi.f90)
target_link_libraries(mpi_f_stub PRIVATE elapse)
target_include_directories(mpi_f_stub INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(mpi_f_stub PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
add_library(MPI::MPI_Fortran ALIAS mpi_f_stub)

# add_library(mpi_f_legacy mpi_legacy.f90)  # use MUMPS::mpiseq to avoid symbol clashes

add_library(scalapack scalapack.f90 blacs.f90)
add_library(SCALAPACK::SCALAPACK ALIAS scalapack)
