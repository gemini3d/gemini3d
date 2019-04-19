# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindBLACS
---------

Finds the BLACSS library


Result Variables
^^^^^^^^^^^^^^^^

BLACS_LIBRARIES
  libraries to be linked

BLACS_INCLUDE_DIRS
  dirs to be included

 Note: choosing static compilation of Scalapack means that BLACS can be entirely within libscalapack.a

#]=======================================================================]

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
            NAMES blacsF77init blacsF77init-openmpi)

find_library(BLACS_LIBRARY
            NAMES blacs blacs-mpi blacs-openmpi)

if(NOT BLACS_LIBRARY)
  return()
endif()

if(BLACS_CINIT)
  list(APPEND BLACS_LIBRARY ${BLACS_CINIT})
endif()

if(BLACS_INIT)
  list(APPEND BLACS_LIBRARY ${BLACS_INIT})
endif()

end()

endfunction(getlibs)


if(NOT DEFINED BLACS_FIND_COMPONENTS)
  set(BLACS_FIND_COMPONENTS OpenMPI)
endif()

getlibs()

if(BLACS_LIBRARY)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${BLACS_LIBRARY} MPI::MPI_Fortran)
  check_fortran_function_exists(blacs_gridmap BLACS_OK)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLACS
  REQUIRED_VARS BLACS_LIBRARY BLACS_OK
  HANDLE_COMPONENTS)


if(BLACS_FOUND)
  set(BLACS_LIBRARIES ${BLACS_LIBRARY})
endif()

mark_as_advanced(BLACS_LIBRARY)
