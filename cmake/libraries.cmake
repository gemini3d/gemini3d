# -------  Libraries-----------
# MKL
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
find_package(MKL COMPONENTS MPI)
set(MKLROOT $ENV{MKLROOT})
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  include_directories(${MKL_INCLUDE_DIRS})
  set(BLA_VENDOR Intel)
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# MPI
find_package(MPI REQUIRED)
add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})

message(STATUS "MPI include: " ${MPI_Fortran_INCLUDE_DIRS})
message(STATUS "MPI lib: " ${MPI_Fortran_LIBRARIES})
message(STATUS "MPI flags: " ${MPI_Fortran_COMPILE_OPTIONS})
message("MPI compiler: " ${MPI_Fortran_COMPILER})


# SCALAPACK
list(APPEND SCALAPACK_ROOT ${LIB_DIR}/scalapack)
find_package(SCALAPACK REQUIRED)

# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
list(APPEND METIS_ROOT ${LIB_DIR}/metis)
list(APPEND Scotch_ROOT ${LIB_DIR}/scotch)
list(APPEND MUMPS_ROOT ${LIB_DIR}/MUMPS)

# "d" is float64  "s" is float32
find_package(MUMPS REQUIRED COMPONENTS d)  # apt install libmumps-dev
