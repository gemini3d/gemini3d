# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindHWLOC
-------
Michael Hirsch, Ph.D.

Finds the hwloc library, required by OpenMPI and also useful by itself.

Imported Targets
^^^^^^^^^^^^^^^^

HWLOC::HWLOC

Result Variables
^^^^^^^^^^^^^^^^

HWLOC_LIBRARIES
  libraries to be linked

HWLOC_INCLUDE_DIRS
  dirs to be included

#]=======================================================================]

find_package(PkgConfig)
pkg_check_modules(pc_hwloc hwloc)


find_library(HWLOC_LIBRARY
             NAMES hwloc
             HINTS ${pc_hwloc_LIBRARY_DIRS} ${pc_hwloc_LIBDIR})

find_path(HWLOC_INCLUDE_DIR
          NAMES hwloc.h
          HINTS ${pc_hwloc_INCLUDE_DIRS})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC
    REQUIRED_VARS HWLOC_LIBRARY HWLOC_INCLUDE_DIR)

if(HWLOC_FOUND)
# need if _FOUND guard to allow project to autobuild; can't overwrite imported target even if bad
set(HWLOC_LIBRARIES ${HWLOC_LIBRARY})
set(HWLOC_INCLUDE_DIRS ${HWLOC_INCLUDE_DIR})

if(NOT TARGET HWLOC::HWLOC)
  add_library(HWLOC::HWLOC INTERFACE IMPORTED)
  set_target_properties(HWLOC::HWLOC PROPERTIES
                        INTERFACE_LINK_LIBRARIES "${HWLOC_LIBRARY}"
                        INTERFACE_INCLUDE_DIRECTORIES "${HWLOC_INCLUDE_DIR}"
                      )
endif()
endif(HWLOC_FOUND)

mark_as_advanced(HWLOC_INCLUDE_DIR HWLOC_LIBRARY)
