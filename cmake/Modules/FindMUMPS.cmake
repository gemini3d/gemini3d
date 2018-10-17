# - Try to find MUMPS
# https://cmake.org/cmake/help/v3.11/manual/cmake-developer.7.html#find-modules

#.rst:
# FindMUMPS
# -------
#
# Finds the MUMPS library
#
# This will define the following variables::
#
#   MUMPS_FOUND    - True if the system has the MUMPS library
#   MUMPS_VERSION  - The version of the MUMPS library which was found
#   MUMPS_LIBRARIES - the libraries to be linked
#   MUMPS_INCLUDE_DIRS - dirs to be included
#
# Inputs:
#  COMPONENTS  s d c z   list one or more

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  find_package(MKL COMPONENTS MPI)
else()
  find_package(LAPACK)
  if (NOT LAPACK_FOUND)
    return()
  endif()

  # ----- MPI
  find_package(MPI COMPONENTS Fortran)
  if (NOT MPI_FOUND)
    return()
  endif()
endif()

find_package(Threads REQUIRED)
add_compile_options(${MPI_Fortran_COMPILE_OPTIONS})
include_directories(${MPI_Fortran_INCLUDE_DIRS})

#-- METIS
find_package(METIS)
#-- Scotch
find_package(Scotch COMPONENTS ESMUMPS)
#-- Scalapack
find_package(SCALAPACK)
#-------------------------


find_package(PkgConfig)
pkg_check_modules(PC_MUMPS QUIET MUMPS)


find_path(MUMPS_INCLUDE_DIR
          NAMES mumps_compat.h
          PATHS ${PC_MUMPS_INCLUDE_DIRS}
          PATH_SUFFIXES MUMPS include
          HINTS ${MUMPS_ROOT})

find_library(MUMPS_COMMON
             NAMES mumps_common
             PATHS ${PC_MUMPS_LIBRARY_DIRS}
             PATH_SUFFIXES MUMPS lib
             HINTS ${MUMPS_ROOT})

find_library(PORD 
              NAMES pord 
              PATH_SUFFIXES MUMPS lib
              HINTS ${MUMPS_ROOT})
        
unset(MUMPS_LIBRARIES)

FOREACH(comp ${MUMPS_FIND_COMPONENTS})
  find_library(MUMPS_${comp}_lib
              NAMES ${comp}mumps 
              PATH_SUFFIXES MUMPS lib
              HINTS ${MUMPS_ROOT})
  
#  message(STATUS "MUMPS finding " ${comp}mumps " in " MUMPS_${thiscomp})
  if(MUMPS_${comp}_lib)
    list(APPEND MUMPS_LIBRARIES ${MUMPS_${comp}_lib})
    mark_as_advanced(MUMPS_${comp}_lib)
  else()
      message(FATAL_ERROR "did not find" ${MUMPS_${comp}_lib})
  endif()
ENDFOREACH()


set(MUMPS_VERSION ${PC_MUMPS_VERSION})

list(APPEND MUMPS_LIBRARIES ${MUMPS_COMMON} ${PORD})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
    FOUND_VAR MUMPS_FOUND
    REQUIRED_VARS MUMPS_LIBRARIES MUMPS_INCLUDE_DIR
    VERSION_VAR MUMPS_VERSION)

# in this order!
set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR} )

if(METIS_FOUND)
  list(APPEND MUMPS_LIBRARIES ${METIS_LIBRARIES})
endif()
if(Scotch_FOUND)
  list(APPEND MUMPS_LIBRARIES ${Scotch_LIBRARIES})
  list(APPEND MUMPS_INCLUDE_DIRS  ${Scotch_INCLUDE_DIRS})
endif()
if(SCALAPACK_FOUND)
  list(APPEND MUMPS_LIBRARIES ${SCALAPACK_LIBRARIES})
  list(APPEND MUMPS_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIRS})
endif()

list(APPEND MUMPS_LIBRARIES ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) 

list(APPEND MUMPS_INCLUDE_DIRS ${LAPACK_INCLUDE_DIRS}  ${MPI_Fortran_INCLUDE_PATH})

set(MUMPS_DEFINITIONS  ${PC_MUMPS_CFLAGS_OTHER})

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)

