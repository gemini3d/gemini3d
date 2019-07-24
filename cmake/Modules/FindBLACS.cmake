# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindBLACS
---------

Finds the OpenMPI BLACS library.
If you want this via MKL, just use ``find_package(SCALAPACK)``



Result Variables
^^^^^^^^^^^^^^^^

BLACS_LIBRARIES
  libraries to be linked

BLACS_INCLUDE_DIRS
  dirs to be included

 Note: choosing static compilation of Scalapack means that BLACS can be entirely within libscalapack.a

#]=======================================================================]

cmake_policy(VERSION 3.3)

function(getlibs)

if(MPICH IN_LIST BLACS_FIND_COMPONENTS)
  find_library(BLACS_LIBRARY
               NAMES blacs-mpich blacs-mpich2)
elseif(LAM IN_LIST BLACS_FIND_COMPONENTS)
  find_library(BLACS_LIBRARY
               NAMES blacs-lam)
elseif(PVM IN_LIST BLACS_FIND_COMPONENTS)
  find_library(BLACS_LIBRARY
               NAMES blacs-pvm)
elseif(OpenMPI IN_LIST BLACS_FIND_COMPONENTS)

  find_library(BLACS_INIT
    NAMES blacsF77init blacsF77init-openmpi
    PATHS ${SCALAPACK_ROOT})

  # this is the only lib that scalapack/blacs/src provides
  find_library(BLACS_CINIT
    NAMES blacsCinit blacsCinit-openmpi
    PATHS ${SCALAPACK_ROOT})

  find_library(BLACS_LIB
    NAMES blacs blacs-mpi blacs-openmpi
    PATHS ${SCALAPACK_ROOT})

  if(BLACS_LIB)
    list(APPEND BLACS_LIBRARY ${BLACS_LIB})
  endif()

  if(BLACS_CINIT)
    list(APPEND BLACS_LIBRARY ${BLACS_CINIT})
  endif()

  if(BLACS_INIT)
    list(APPEND BLACS_LIBRARY ${BLACS_INIT})
  endif()

  if(BLACS_LIBRARY)
    set(BLACS_OpenMPI_FOUND true PARENT_SCOPE)
  endif()
else()
  message(FATAL_ERROR "Must specify desired COMPONENTS for find_package(BLACS)")
endif()


set(BLACS_LIBRARY ${BLACS_LIBRARY} PARENT_SCOPE)

endfunction(getlibs)

# == main

if(NOT BLACS_FIND_COMPONENTS)
  set(BLACS_FIND_COMPONENTS OpenMPI)
endif()

getlibs()


if(BLACS_LIBRARY)
  include(CheckFortranSourceCompiles)
  set(CMAKE_REQUIRED_LIBRARIES ${BLACS_LIBRARY})
  check_fortran_source_compiles("integer contxt; CALL BLACS_GET(0, 0, CONTXT); end" BLACS_OK SRC_EXT f90)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLACS
  REQUIRED_VARS BLACS_LIBRARY
  HANDLE_COMPONENTS)

if(BLACS_FOUND)
  set(BLACS_LIBRARIES ${BLACS_LIBRARY})
endif()

mark_as_advanced(BLACS_LIBRARY)
