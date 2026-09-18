# ref: https://jcbhmr.me/blog/cosmocc-cmake
# https://gitlab.kitware.com/cmake/cmake/-/issues/25578
#
# Usage:
# cmake --workflow cosmo

set(CMAKE_SYSTEM_NAME Generic)
unset(CMAKE_SYSTEM_PROCESSOR)

# Allow users to provide COSMO_ROOT via environment without editing presets.
if(NOT COSMO_ROOT AND DEFINED ENV{COSMO_ROOT})
  set(COSMO_ROOT "$ENV{COSMO_ROOT}")
endif()

if(DEFINED ENV{CC} AND EXISTS "$ENV{CC}")
  set(CMAKE_C_COMPILER "$ENV{CC}")
else()
  find_program(CMAKE_C_COMPILER
  NAMES cosmocc
  REQUIRED
  PATH_SUFFIXES bin
  HINTS ${COSMO_ROOT}
  )
endif()

if(DEFINED ENV{CXX} AND EXISTS "$ENV{CXX}")
  set(CMAKE_CXX_COMPILER "$ENV{CXX}")
else()
  find_program(CMAKE_CXX_COMPILER
  NAMES cosmoc++
  REQUIRED
  PATH_SUFFIXES bin
  HINTS ${COSMO_ROOT}
  )
endif()

if(DEFINED ENV{ASM} AND EXISTS "$ENV{ASM}")
  set(CMAKE_ASM_COMPILER "$ENV{ASM}")
else()
  find_program(CMAKE_ASM_COMPILER
  NAMES cosmocc
  REQUIRED
  PATH_SUFFIXES bin
  HINTS ${COSMO_ROOT}
  )
endif()


if(DEFINED ENV{AR} AND EXISTS "$ENV{AR}")
  set(CMAKE_AR "$ENV{AR}")
else()
  find_program(CMAKE_AR
  NAMES cosmoar
  REQUIRED
  PATH_SUFFIXES bin
  HINTS ${COSMO_ROOT}
  )
endif()

if(DEFINED ENV{RANLIB} AND EXISTS "$ENV{RANLIB}")
  set(CMAKE_RANLIB "$ENV{RANLIB}")
else()
  find_program(CMAKE_RANLIB
  NAMES cosmoranlib
  REQUIRED
  PATH_SUFFIXES bin
  HINTS ${COSMO_ROOT}
  )
endif()

set(CMAKE_CXX_FLAGS_INIT "-fexceptions -frtti")

set(COSMOPOLITAN 1)
set(UNIX 1)

# Cosmopolitan uses only static libraries and no RPATH
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_SKIP_RPATH ON)

# The applications can run on the host platform
set(CMAKE_CROSSCOMPILING OFF)

set(ffilesystem_fortran false)
