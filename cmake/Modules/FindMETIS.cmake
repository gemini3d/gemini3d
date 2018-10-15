# - Try to find METIS (not ParMETIS) library
# https://cmake.org/cmake/help/v3.11/manual/cmake-developer.7.html#find-modules

#.rst:
# FindMETIS
# -------
# Michael Hirsch, Ph.D.
#
# Finds the METIS library
#
# This will define the following variables::
#
#   METIS_FOUND    - True if the system has the METIS library
#   METIS_VERSION  - The version of the METIS library which was found

find_package(PkgConfig)
pkg_check_modules(PC_METIS QUIET METIS)


find_library(METIS_LIBRARY
             NAMES metis
             PATHS ${PC_METIS_LIBRARY_DIRS}
             PATH_SUFFIXES METIS lib libmetis build/Linux-x86_64/libmetis
             HINTS ${METIS_ROOT})

find_path(METIS_INCLUDE_DIR
          NAMES metis.h
          PATHS ${PC_METIS_INCLUDE_DIRS}
          PATH_SUFFIXES METIS include
          HINTS ${METIS_ROOT})

set(METIS_VERSION ${PC_METIS_VERSION})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
    FOUND_VAR METIS_FOUND
    REQUIRED_VARS METIS_LIBRARY METIS_INCLUDE_DIR
    VERSION_VAR METIS_VERSION)

if(METIS_FOUND)
  set(METIS_LIBRARIES ${METIS_LIBRARY})
  set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
  set(METIS_DEFINITIONS  ${PC_METIS_CFLAGS_OTHER})
endif()


mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY)

