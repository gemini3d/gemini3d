add_library(mpi_c_stub mpic.c)
target_include_directories(mpi_c_stub INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(elapse OBJECT elapse.c)
target_include_directories(elapse INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(mpi_f_stub mpi.f90 $<TARGET_OBJECTS:elapse>)

# --- what you'll link to

add_library(MPI::MPI_C INTERFACE IMPORTED GLOBAL)
target_link_libraries(MPI::MPI_C INTERFACE mpi_c_stub)

add_library(MPI::MPI_Fortran INTERFACE IMPORTED GLOBAL)
target_link_libraries(MPI::MPI_Fortran INTERFACE mpi_f_stub)

# add_library(mpi_f_legacy mpi_legacy.f90)  # use MUMPS::MPISEQ to avoid symbol clashes

add_library(scalapack scalapack.f90 blacs.f90)
