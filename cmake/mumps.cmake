# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

unset(_mumps_extra)

# BLACS
if(BLACS_ROOT)
  find_package(BLACS REQUIRED)
  list(APPEND _mumps_extra ${BLACS_LIBRARIES})
endif()

# Metis
if(METIS_ROOT)
  find_package(METIS REQUIRED)
  list(APPEND _mumps_extra ${METIS_LIBRARIES})
endif()

# Scotch
if(Scotch_ROOT)
  find_package(Scotch REQUIRED COMPONENTS ESMUMPS)
  list(APPEND _mumps_extra ${Scotch_LIBRARIES})
endif()

# SCALAPACK always needed for MUMPS
if(DEFINED ENV{MKLROOT})
  find_package(SCALAPACK REQUIRED COMPONENTS MKL)
else()
  find_package(SCALAPACK REQUIRED COMPONENTS OpenMPI)
endif()

# LAPACK
if(DEFINED ENV{MKLROOT})
  find_package(LAPACK REQUIRED COMPONENTS MKL)
else()
  find_package(LAPACK REQUIRED)
endif()

# -- verify Scalapack links

set(CMAKE_REQUIRED_INCLUDES ${SCALAPACK_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})

include(CheckFortranSourceRuns)
file(READ ${CMAKE_SOURCE_DIR}/tests/test_scalapack.f90 _code)
check_fortran_source_runs(${_code} SCALAPACK_OK SRC_EXT f90)

if(NOT SCALAPACK_OK)
message(FATAL_ERROR "Scalapack ${SCALAPACK_LIBRARIES} not working with LAPACK ${LAPACK_LIBRARIES} and ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()


# --- MUMPS
if(realbits EQUAL 64)
  set(_mumpscomp d)
elseif(realbits EQUAL 32)
  set(_mumpscomp s)
else()
  message(FATAL_ERROR "MUMPS has only real32, real64")
endif()

find_package(MUMPS REQUIRED COMPONENTS ${_mumpscomp})

# rather than appending this libraries everywhere, just put them together here.
list(APPEND MUMPS_LIBRARIES ${SCALAPACK_LIBRARIES})
list(APPEND MUMPS_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIRS})

list(APPEND MUMPS_LIBRARIES ${_mumps_extra})

if(MUMPS_ROOT)  # not a system library, need lapack
  list(APPEND MUMPS_LIBRARIES ${LAPACK_LIBRARIES})
endif()

# -- verify MUMPS works

include(CheckFortranSourceCompiles)

set(CMAKE_REQUIRED_LIBRARIES ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
set(CMAKE_REQUIRED_INCLUDES ${MUMPS_INCLUDE_DIRS})

check_fortran_source_compiles("include '${_mumpscomp}mumps_struc.h'
type(${_mumpscomp}mumps_struc) :: mumps_par
end"
  MUMPS_OK SRC_EXT f90)

if(NOT MUMPS_OK)
message(FATAL_ERROR "MUMPS ${MUMPS_LIBRARIES} not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
endif()
