# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:

FindSCALAPACK
-------------

authored by SciVision: www.scivision.dev

Finds SCALAPACK libraries for MKL, OpenMPI and MPICH.
Intel MKL relies on having environment variable MKLROOT set, typically by sourcing
mklvars.sh beforehand.

This module does NOT find LAPACK.

COMPONENTS
^^^^^^^^^^

``MKL``
  Intel MKL for MSVC, oneAPI, GCC.
  Working with IntelMPI (default Window, Linux), MPICH (default Mac) or OpenMPI (Linux only).
``MKL64``
  MKL 64-bit integers  (default is 32-bit integers)
``TBB``
  MKL only: Intel MPI + TBB (default is sequential)
``OpenMP``
  MKL only: use OpenMP (default is sequential)

``AOCL``
  AMD ScaLAPACK fork of Netlib ScaLAPACK.
  Requires LAPACK AOCL
  https://www.amd.com/en/developer/aocl/scalapack.html
``AOCL64``
  AOCL 64-bit integers  (default is 32-bit integers)

``STATIC``
  Library search default on non-Windows is shared then static. On Windows default search is static only.
  Specifying STATIC component searches for static libraries only.

Result Variables
^^^^^^^^^^^^^^^^

``SCALAPACK_FOUND``
  SCALapack libraries were found
``SCALAPACK_<component>_FOUND``
  SCALAPACK <component> specified was found
``SCALAPACK_LIBRARIES``
  SCALapack library files
``SCALAPACK_INCLUDE_DIRS``
  SCALapack include directories


References
^^^^^^^^^^

* Pkg-Config and MKL:  https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-and-pkg-config-tool
* MKL for Windows: https://software.intel.com/en-us/mkl-windows-developer-guide-static-libraries-in-the-lib-intel64-win-directory
* MKL Windows directories: https://software.intel.com/en-us/mkl-windows-developer-guide-high-level-directory-structure
* MKL link-line advisor: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#]=======================================================================]

include(CheckSourceCompiles)

set(SCALAPACK_LIBRARY)  # avoids appending to prior FindScalapack

#===== functions

function(scalapack_check)

# some OpenMPI builds need -pthread
find_package(Threads)


set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LINK_OPTIONS)
set(CMAKE_REQUIRED_INCLUDES ${SCALAPACK_INCLUDE_DIR} ${LAPACK_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARY})
if(BLACS_LIBRARY)
  list(APPEND CMAKE_REQUIRED_LIBRARIES ${BLACS_LIBRARY})
endif()
list(APPEND CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(STATIC IN_LIST SCALAPACK_FIND_COMPONENTS AND
  NOT WIN32 AND
  MKL IN_LIST SCALAPACK_FIND_COMPONENTS AND
  CMAKE_VERSION VERSION_GREATER_EQUAL 3.24
  )
  set(CMAKE_REQUIRED_LIBRARIES $<LINK_GROUP:RESCAN,${CMAKE_REQUIRED_LIBRARIES}>)
endif()
# MPI needed for IntelLLVM

check_source_compiles(Fortran
"program test
use, intrinsic :: iso_fortran_env, only : real64
implicit none
real(real64), external :: pdlamch
integer :: ictxt
print *, pdlamch(ictxt, 'E')
end program"
SCALAPACK_d_FOUND
)

check_source_compiles(Fortran
"program test
use, intrinsic :: iso_fortran_env, only : real32
implicit none
real(real32), external :: pslamch
integer :: ictxt
print *, pslamch(ictxt, 'E')
end program"
SCALAPACK_s_FOUND
)

if(SCALAPACK_s_FOUND OR SCALAPACK_d_FOUND)
  set(SCALAPACK_links true PARENT_SCOPE)
endif()

endfunction()


macro(scalapack_mkl)
# https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2025-2/cmake-config-for-onemkl.html

set(ENABLE_SCALAPACK true)
set(ENABLE_BLAS true)

set(MKL_INTERFACE "lp64")
if(MKL64 IN_LIST SCALAPACK_FIND_COMPONENTS)
  string(PREPEND MKL_INTERFACE "i")
endif()

# MKL_THREADING default: "intel_thread" which is Intel OpenMP
# some systems have messed up OpenMP, so sequential unless requested
if(NOT DEFINED MKL_THREADING)
  if(TBB IN_LIST SCALAPACK_FIND_COMPONENTS)
    set(MKL_THREADING "tbb_thread")
  elseif(OpenMP IN_LIST SCALAPACK_FIND_COMPONENTS)
    set(MKL_THREADING "intel_thread")
  else()
    set(MKL_THREADING "sequential")
  endif()
endif()

set(MKL_SYCL_MPI false)
set(MKL_SYCL_LINK false)
# for Intel oneAPI 2025.2, we don't need SYCL

# default: dynamic
if(STATIC IN_LIST SCALAPACK_FIND_COMPONENTS)
  set(MKL_LINK "static")
endif()

find_package(MKL CONFIG)

if(NOT MKL_FOUND)
  return()
endif()

# get_property(SCALAPACK_COMPILE_OPTIONS TARGET MKL::MKL PROPERTY INTERFACE_COMPILE_OPTIONS)
# flags are empty generator expressions that trip up check_source_compiles

get_property(SCALAPACK_INCLUDE_DIR TARGET MKL::MKL PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
get_property(SCALAPACK_LIBRARY TARGET MKL::MKL PROPERTY INTERFACE_LINK_LIBRARIES)

set(SCALAPACK_MKL_FOUND true)

foreach(c IN ITEMS TBB MKL64 OpenMP)
  if(${c} IN_LIST SCALAPACK_FIND_COMPONENTS)
    set(SCALAPACK_${c}_FOUND true)
  endif()
endforeach()

endmacro()

#==========================

function(scalapack_aocl)

set(_nodef_scalapack)
if(DEFINED SCALAPACK_ROOT)
  set(_nodef_scalapack NO_DEFAULT_PATH)
endif()

set(_s "LP64")
if(AOCL64 IN_LIST SCALAPACK_FIND_COMPONENTS)
  string(PREPEND _s "I")
endif()

find_library(SCALAPACK_LIBRARY
NAMES scalapack
PATH_SUFFIXES lib/${_s}
HINTS ${SCALAPACK_ROOT} $ENV{SCALAPACK_ROOT}
${_nodef_scalapack}
DOC "AOCL SCALAPACK library"
)

if(SCALAPACK_LIBRARY)
  set(SCALAPACK_AOCL_FOUND true PARENT_SCOPE)
endif()

endfunction()

#===========================

function(scalapack_netlib)

if(BUILD_SHARED_LIBS)
  set(_s shared)
else()
  set(_s static)
endif()
list(APPEND _s openmpi/lib mpich/lib)

find_library(SCALAPACK_LIBRARY
NAMES scalapack scalapack-openmpi scalapack-mpich
NAMES_PER_DIR
PATH_SUFFIXES ${_s}
DOC "SCALAPACK library"
)

# some systems have libblacs as a separate file, instead of being subsumed in libscalapack.
if(NOT DEFINED BLACS_ROOT)
  cmake_path(GET SCALAPACK_LIBRARY PARENT_PATH BLACS_ROOT)
endif()

find_library(BLACS_LIBRARY
NAMES blacs
NO_DEFAULT_PATH
HINTS ${BLACS_ROOT}
DOC "BLACS library"
)

endfunction()

# === main

if(NOT DEFINED SCALAPACK_CRAY AND DEFINED ENV{CRAYPE_VERSION})
  set(SCALAPACK_CRAY true)
endif()

if(NOT SCALAPACK_CRAY)
  if(NOT MKL IN_LIST SCALAPACK_FIND_COMPONENTS AND DEFINED ENV{MKLROOT} AND IS_DIRECTORY "$ENV{MKLROOT}")
    list(APPEND SCALAPACK_FIND_COMPONENTS MKL)
  endif()
endif()

if(STATIC IN_LIST SCALAPACK_FIND_COMPONENTS)
  set(_orig_suff ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

if(MKL IN_LIST SCALAPACK_FIND_COMPONENTS OR MKL64 IN_LIST SCALAPACK_FIND_COMPONENTS)
  scalapack_mkl()
elseif(SCALAPACK_CRAY)
  # Cray PE has Scalapack build into LibSci. Use Cray compiler wrapper.
elseif(AOCL IN_LIST LAPACK_FIND_COMPONENTS)
  scalapack_aocl()
else()
  scalapack_netlib()
endif()

if(STATIC IN_LIST SCALAPACK_FIND_COMPONENTS)
  if(SCALAPACK_LIBRARY)
    set(SCALAPACK_STATIC_FOUND true)
  endif()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_orig_suff})
endif()

# --- Check that Scalapack links

if(SCALAPACK_CRAY OR SCALAPACK_LIBRARY)
  scalapack_check()
endif()

# --- Finalize

include(FindPackageHandleStandardArgs)

if(SCALAPACK_CRAY)
  find_package_handle_standard_args(SCALAPACK HANDLE_COMPONENTS
  REQUIRED_VARS SCALAPACK_links
  )
else()
  find_package_handle_standard_args(SCALAPACK HANDLE_COMPONENTS
  REQUIRED_VARS SCALAPACK_LIBRARY SCALAPACK_links
  )
endif()

if(SCALAPACK_FOUND)
  # need if _FOUND guard as can't overwrite imported target even if bad
  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
  if(BLACS_LIBRARY)
    list(APPEND SCALAPACK_LIBRARIES ${BLACS_LIBRARY})
  endif()

  set(SCALAPACK_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIR})

  message(VERBOSE "Scalapack libraries: ${SCALAPACK_LIBRARIES}
Scalapack include directories: ${SCALAPACK_INCLUDE_DIRS}")

  if(NOT TARGET SCALAPACK::SCALAPACK)
    add_library(SCALAPACK::SCALAPACK INTERFACE IMPORTED)
    set_property(TARGET SCALAPACK::SCALAPACK PROPERTY INTERFACE_LINK_LIBRARIES "${SCALAPACK_LIBRARIES}")
    set_property(TARGET SCALAPACK::SCALAPACK PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${SCALAPACK_INCLUDE_DIR}")

    # For MKL, we don't use FindLapack, so define LAPACK::LAPACK as alias
    if(MKL_FOUND AND NOT TARGET LAPACK::LAPACK)
      add_library(LAPACK::LAPACK ALIAS SCALAPACK::SCALAPACK)
    endif()
  endif()
endif()

mark_as_advanced(SCALAPACK_LIBRARY SCALAPACK_INCLUDE_DIR)
