# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindBLACS
---------

Finds the BLACSS library



Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

BLACS::BLACS
  the BLASC library components requested.

Result Variables
^^^^^^^^^^^^^^^^

BLACS_LIBRARIES
  libraries to be linked

BLACS_INCLUDE_DIRS
  dirs to be included

 Note: choosing static compilation of Scalapack means that BLACS can be entirely within libscalapack.a

#]=======================================================================]

find_library(BLACS_LIBRARY
            NAMES blacs blacs-pvm blacs-mpi blacs-openmpi blacsF77init-openmpi blacs-mpich blacs-mpich2 blacs-lam
            PATH_SUFFIXES lib)


find_library(BLACS_OPENMPI 
            NAMES blacs-openmpi
            PATH_SUFFIXES lib)

find_library(BLACS_CINIT 
            NAMES blacsCinit-openmpi
            PATH_SUFFIXES lib)

list(APPEND BLACS_LIBRARY ${BLACS_OPENMPI} ${BLACS_CINIT})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLACS
    REQUIRED_VARS BLACS_LIBRARY)  
# don't put BLACS_LIBRARY REQUIRED_VARS because it might be in libBLACS.a)

if(BLACS_FOUND)
  set(BLACS_LIBRARIES ${BLACS_LIBRARY})
endif()

mark_as_advanced(BLACS_LIBRARY)
