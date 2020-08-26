add_library(mpi_c_stub mpic.c)
target_include_directories(mpi_c_stub INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_library(MPI::MPI_C ALIAS mpi_c_stub)

add_library(elapse OBJECT elapse.c)

add_library(mpi_f_stub mpi.f90)
target_link_libraries(mpi_f_stub PRIVATE elapse)
add_library(MPI::MPI_Fortran ALIAS mpi_f_stub)

# add_library(mpi_f_legacy mpi_legacy.f90)  # use mumps::mpiseq to avoid symbol clashes

add_library(scalapack scalapack.f90 blacs.f90)
add_library(SCALAPACK::SCALAPACK ALIAS scalapack)

foreach(t mpi_f_stub mpi_c_stub scalapack)
  set_target_properties(${t} PROPERTIES
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach()
