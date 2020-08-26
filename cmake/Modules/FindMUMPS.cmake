# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindMUMPS
---------

Finds the MUMPS library.
Note that MUMPS generally requires SCALAPACK and LAPACK as well.
PORD is always used, in addition to the optional METIS or SCOTCH, which would be found externally.

COMPONENTS
  s d c z   list one or more. Defaults to ``d``.
  mpiseq    for -Dparallel=false, a stub MPI & Scalapack library

Result Variables
^^^^^^^^^^^^^^^^

MUMPS_LIBRARIES
  libraries to be linked

MUMPS_INCLUDE_DIRS
  dirs to be included

#]=======================================================================]

set(MUMPS_LIBRARY)  # don't endlessly append

function(mumps_libs)

find_path(MUMPS_INCLUDE_DIR
          NAMES mumps_compat.h
          DOC "MUMPS common header")
if(NOT MUMPS_INCLUDE_DIR)
  return()
endif()

find_library(MUMPS_COMMON
             NAMES mumps_common
             NAMES_PER_DIR
             DOC "MUMPS common libraries")
if(NOT MUMPS_COMMON)
  return()
endif()

find_library(PORD
             NAMES pord
             NAMES_PER_DIR
             DOC "simplest MUMPS ordering library")
if(NOT PORD)
  return()
endif()

if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  find_library(MUMPS_mpiseq_LIB
               NAMES mpiseq
               NAMES_PER_DIR
               DOC "No-MPI stub library")
  if(NOT MUMPS_mpiseq_LIB)
    message(WARNING "MUMPS mpiseq not found")
    return()
  endif()

  find_path(MUMPS_mpiseq_INC
            NAMES mpif.h
            DOC "MUMPS mpiseq header")
  if(NOT MUMPS_mpiseq_INC)
    message(WARNING "MUMPS mpiseq mpif.h not found")
    return()
  endif()

  set(MUMPS_mpiseq_FOUND true PARENT_SCOPE)
  set(MUMPS_mpiseq_LIB ${MUMPS_mpiseq_LIB} PARENT_SCOPE)
  set(MUMPS_mpiseq_INC ${MUMPS_mpiseq_INC} PARENT_SCOPE)
endif()

foreach(comp ${MUMPS_FIND_COMPONENTS})
  if(comp STREQUAL mpiseq)
    continue()
  endif()

  find_library(MUMPS_${comp}_lib
               NAMES ${comp}mumps
               NAMES_PER_DIR)

  if(NOT MUMPS_${comp}_lib)
    message(WARNING "MUMPS ${comp} not found")
    return()
  endif()

  set(MUMPS_${comp}_FOUND true PARENT_SCOPE)
  list(APPEND MUMPS_LIBRARY ${MUMPS_${comp}_lib})
endforeach()

set(MUMPS_LIBRARY ${MUMPS_LIBRARY} ${MUMPS_COMMON} ${PORD} PARENT_SCOPE)
set(MUMPS_INCLUDE_DIR ${MUMPS_INCLUDE_DIR} PARENT_SCOPE)

endfunction(mumps_libs)


if(NOT MUMPS_FIND_COMPONENTS)
  list(APPEND MUMPS_FIND_COMPONENTS d)
endif()

mumps_libs()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS
  REQUIRED_VARS MUMPS_LIBRARY MUMPS_INCLUDE_DIR
  HANDLE_COMPONENTS)

if(MUMPS_FOUND)
# need if _FOUND guard to allow project to autobuild; can't overwrite imported target even if bad
set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  list(APPEND MUMPS_LIBRARIES ${MUMPS_mpiseq_LIB})
  list(APPEND MUMPS_INCLUDE_DIRS ${MUMPS_mpiseq_INC})
endif()

if(NOT TARGET MUMPS::MUMPS)
  add_library(MUMPS::MUMPS INTERFACE IMPORTED)
  set_target_properties(MUMPS::MUMPS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${MUMPS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_INCLUDE_DIR}")
endif()

if(mpiseq IN_LIST MUMPS_FIND_COMPONENTS)
  if(NOT TARGET MUMPS::mpiseq)
    add_library(MUMPS::mpiseq INTERFACE IMPORTED)
    set_target_properties(MUMPS::mpiseq PROPERTIES
      INTERFACE_LINK_LIBRARIES "${MUMPS_mpiseq_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUMPS_mpiseq_INC}")
  endif()
endif()

endif(MUMPS_FOUND)

mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)
