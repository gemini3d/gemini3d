# - Try to find Netlib LAPACK95 library
# https://cmake.org/cmake/help/v3.11/manual/cmake-developer.7.html#find-modules

#.rst:
# FindLAPACK95
# -------
# Michael Hirsch, Ph.D.
#
# Finds the LAPACK95 static library  (need the .a/.lib specified or CMake won't find it, unlike .so/.dll that can be omitted)
# tested with Netlib LAPACK95
#
# This will define the following variables::
#
#   LAPACK95_FOUND    - True if the system has the LAPACK95 library
#   LAPACK95_VERSION  - The version of the LAPACK95 library which was found

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

  find_path(LAPACK95_INCLUDE_DIR 
            NAMES lapack95.mod
            HINTS $ENV{MKLROOT}/include 
                  $ENV{MKLROOT}/include/intel64/lp64)
          
  foreach(slib mkl_blas95_lp64 mkl_lapack95_lp64 mkl_intel_lp64 mkl_sequential mkl_core)
    find_library(LAPACK95_${slib}_LIBRARY
             NAMES ${slib}
             HINTS $ENV{MKLROOT}/lib
                   $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)
    if(NOT LAPACK95_${slib}_LIBRARY)
      message(FATAL_ERROR "NOT FOUND: " ${slib} ${LAPACK95_${slib}_LIBRARY})
    endif()
#    message(STATUS "Intel MKL LAPACK95 FOUND: " ${slib} ${LAPACK95_${slib}_LIBRARY})
    list(APPEND LAPACK95_LIBRARY ${LAPACK95_${slib}_LIBRARY})
    mark_as_advanced(LAPACK95_${slib}_LIBRARY)
  endforeach()
  list(APPEND LAPACK95_LIBRARY pthread dl m)

else()
  find_package(LAPACK REQUIRED)


  find_package(PkgConfig)
  pkg_check_modules(PC_LAPACK95 QUIET LAPACK95)


  find_library(LAPACK95_LIBRARY
               NAMES lapack95.a
               PATHS ${PC_LAPACK95_LIBRARY_DIRS}
               PATH_SUFFIXES LAPACK95
               HINTS ${LAPACK95_ROOT})

  find_path(LAPACK95_INCLUDE_DIR
            NAMES f95_lapack.mod
            PATHS ${PC_LAPACK95_INCLUDE_DIRS}
            PATH_SUFFIXES LAPACK95 lapack95_modules
            HINTS ${LAPACK95_ROOT})

  set(LAPACK95_VERSION ${PC_LAPACK95_VERSION})

endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK95
    FOUND_VAR LAPACK95_FOUND
    REQUIRED_VARS LAPACK95_LIBRARY LAPACK95_INCLUDE_DIR
    VERSION_VAR LAPACK95_VERSION)

if(LAPACK95_FOUND)
  set(LAPACK95_LIBRARIES ${LAPACK95_LIBRARY} ${LAPACK_LIBRARIES})  # MUST be in this order!
  set(LAPACK95_INCLUDE_DIRS ${LAPACK95_INCLUDE_DIR} ${LAPACK_INCLUDE_DIRS})
  set(LAPACK95_DEFINITIONS  ${PC_LAPACK95_CFLAGS_OTHER})
endif()

mark_as_advanced(LAPACK95_INCLUDE_DIR LAPACK95_LIBRARY)

