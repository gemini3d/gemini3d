# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:

FindSCALAPACK
-------------

by Michael Hirsch, Ph.D. www.scivision.dev

Finds SCALAPACK libraries for MKL, OpenMPI and MPICH.
Intel MKL relies on having environment variable MKLROOT set, typically by sourcing
mklvars.sh beforehand.
Intended to work with Intel MKL at least through version 2021.

This module does NOT find LAPACK.

Parameters
^^^^^^^^^^

``MKL``
  Intel MKL for MSVC, ICL, ICC, GCC and PGCC. Working with IntelMPI (default Window, Linux), MPICH (default Mac) or OpenMPI (Linux only).

``OpenMPI``
  OpenMPI interface

``MPICH``
  MPICH interface


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

#===== functions

function(mkl_scala)

set(_mkl_libs ${ARGV})

foreach(s ${_mkl_libs})
  find_library(SCALAPACK_${s}_LIBRARY
           NAMES ${s}
           PATHS ENV MKLROOT ENV I_MPI_ROOT ENV TBBROOT
           PATH_SUFFIXES
             lib/intel64 lib/intel64_win
             intel64/lib/release
             lib/intel64/gcc4.7 ../tbb/lib/intel64/gcc4.7
             lib/intel64/vc_mt ../tbb/lib/intel64/vc_mt
             ../compiler/lib/intel64
           HINTS ${MKL_LIBRARY_DIRS} ${MKL_LIBDIR}
           NO_DEFAULT_PATH)
  if(NOT SCALAPACK_${s}_LIBRARY)
    message(WARNING "MKL component not found: " ${s})
    return()
  endif()

  list(APPEND SCALAPACK_LIBRARY ${SCALAPACK_${s}_LIBRARY})
endforeach()


find_path(SCALAPACK_INCLUDE_DIR
  NAMES mkl_scalapack.h
  PATHS ENV MKLROOT ENV I_MPI_ROOT ENV TBBROOT
  PATH_SUFFIXES
    include
    include/intel64/lp64
  HINTS ${MKL_INCLUDE_DIRS})

if(NOT SCALAPACK_INCLUDE_DIR)
  message(WARNING "MKL Include Dir not found")
  return()
endif()

list(APPEND SCALAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIRS})

set(SCALAPACK_MKL_FOUND true PARENT_SCOPE)
set(SCALAPACK_LIBRARY ${SCALAPACK_LIBRARY} PARENT_SCOPE)
set(SCALAPACK_INCLUDE_DIR ${SCALAPACK_INCLUDE_DIR} PARENT_SCOPE)

endfunction(mkl_scala)

# === main

if(NOT (OpenMPI IN_LIST SCALAPACK_FIND_COMPONENTS
        OR MPICH IN_LIST SCALAPACK_FIND_COMPONENTS
        OR MKL IN_LIST SCALAPACK_FIND_COMPONENTS))
if(DEFINED ENV{MKLROOT})
  list(APPEND SCALAPACK_FIND_COMPONENTS MKL)
  if(APPLE)
    list(APPEND SCALAPACK_FIND_COMPONENTS MPICH)
  endif()
else()
  list(APPEND SCALAPACK_FIND_COMPONENTS OpenMPI)
endif()
endif()

find_package(PkgConfig)

# some systems (Ubuntu 16.04) need BLACS explicitly, when it isn't statically compiled into libscalapack
# other systems (homebrew, Ubuntu 18.04) link BLACS into libscalapack, and don't need BLACS as a separately linked library.
if(NOT MKL IN_LIST SCALAPACK_FIND_COMPONENTS)
  find_package(BLACS COMPONENTS ${SCALAPACK_FIND_COMPONENTS})
endif()

if(MKL IN_LIST SCALAPACK_FIND_COMPONENTS)

if(BUILD_SHARED_LIBS)
  set(_mkltype dynamic)
else()
  set(_mkltype static)
endif()

if(WIN32)
  set(_impi impi)
else()
  unset(_impi)
endif()

pkg_check_modules(MKL mkl-${_mkltype}-lp64-iomp)

if(OpenMPI IN_LIST SCALAPACK_FIND_COMPONENTS)
  mkl_scala(mkl_scalapack_lp64 mkl_blacs_openmpi_lp64)
  set(SCALAPACK_OpenMPI_FOUND ${SCALAPACK_MKL_FOUND})
elseif(MPICH IN_LIST SCALAPACK_FIND_COMPONENTS)
  if(APPLE)
    mkl_scala(mkl_scalapack_lp64 mkl_blacs_mpich_lp64)
  elseif(WIN32)
    mkl_scala(mkl_scalapack_lp64 mkl_blacs_mpich2_lp64.lib mpi.lib fmpich2.lib)
  else()  # MPICH linux is just like IntelMPI
    mkl_scala(mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64)
  endif()
  set(SCALAPACK_MPICH_FOUND ${SCALAPACK_MKL_FOUND})
else()
  mkl_scala(mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 ${_impi})
endif()

elseif(OpenMPI IN_LIST SCALAPACK_FIND_COMPONENTS)

pkg_check_modules(SCALAPACK scalapack-openmpi)

find_library(SCALAPACK_LIBRARY
              NAMES scalapack scalapack-openmpi
              HINTS ${SCALAPACK_LIBRARY_DIRS} ${SCALAPACK_LIBDIR})

if(SCALAPACK_LIBRARY)
  set(SCALAPACK_OpenMPI_FOUND true)
endif()

elseif(MPICH IN_LIST SCALAPACK_FIND_COMPONENTS)

pkg_check_modules(SCALAPACK scalapack-mpich)

find_library(SCALAPACK_LIBRARY
              NAMES scalapack-mpich scalapack-mpich2
              HINTS ${SCALAPACK_LIBRARY_DIRS} ${SCALAPACK_LIBDIR})

if(SCALAPACK_LIBRARY)
  set(SCALAPACK_MPICH_FOUND true)
endif()

endif()

# --- Finalize

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK
  REQUIRED_VARS SCALAPACK_LIBRARY
  HANDLE_COMPONENTS)

if(SCALAPACK_FOUND)
  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
  set(SCALAPACK_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIR})

  if(BLACS_FOUND)
    list(APPEND SCALAPACK_LIBRARIES ${BLACS_LIBRARIES})
  endif()
endif()

mark_as_advanced(SCALAPACK_LIBRARY SCALAPACK_INCLUDE_DIR)
