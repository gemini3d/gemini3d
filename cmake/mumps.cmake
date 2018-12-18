# MUMPS -- use cmake -DMUMPS_ROOT= for hint.
#
# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

# Metis
if(METIS_ROOT)
  find_package(METIS)
endif()

# Scotch
if(Scotch_ROOT)
  find_package(Scotch COMPONENTS ESMUMPS)
endif()

# BLACS
if(BLACS_ROOT)
  find_package(BLACS)
endif()

# SCALAPACK always needed for MUMPS
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  find_package(SCALAPACK REQUIRED COMPONENTS IntelSeq)
elseif(DEFINED ENV{MKLROOT})
  find_package(SCALAPACK REQUIRED COMPONENTS MKL)
else()
  find_package(SCALAPACK REQUIRED COMPONENTS OpenMPI)
endif()


# Mumps
if(realbits EQUAL 64)
  set(mumpscomp d)
elseif(realbits EQUAL 32)
  set(mumpscomp s)
else()
  message(FATAL_ERROR "MUMPS has only real32, real64")
endif()

find_package(LAPACK REQUIRED)


find_package(MUMPS REQUIRED COMPONENTS ${mumpscomp})
list(APPEND MUMPS_LIBRARIES ${SCALAPACK_LIBRARIES})

if(BLACS_FOUND)
  list(APPEND MUMPS_LIBRARIES ${BLACS_LIBRARIES})
endif()
  
if(Scotch_FOUND)  
  list(APPEND MUMPS_LIBRARIES ${Scotch_LIBRARIES})
endif()

if(METIS_FOUND)
  list(APPEND MUMPS_LIBRARIES ${METIS_LIBRARIES})
endif()

if(MUMPS_ROOT)
  list(APPEND MUMPS_LIBRARIES ${LAPACK_LIBRARIES})
endif()
