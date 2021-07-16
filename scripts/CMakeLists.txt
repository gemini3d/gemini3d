# for systems where you can't use "cmake -P scripts/install_prereq.cmake"
# this project builds the prereqs

cmake_minimum_required(VERSION 3.20...3.21)

project(BuildGeminiPrereq LANGUAGES C Fortran)

set(arith "s;d")
# NOTE: lists must be passed via CMAKE_CACHE_ARGS using -Dvar:STRING=${var}
# the variable type must be included e.g. STRING

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY EP_UPDATE_DISCONNECTED true)

include(ExternalProject)

include(${PROJECT_SOURCE_DIR}/../cmake/config/libraries.cmake)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(prefix ~/lib_gcc)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(prefix ~/lib_intel)
  else()
    set(prefix ~/lib)
  endif()

  if(CMAKE_VERSION VERSION_LESS 3.21)
    get_filename_component(prefix ${prefix} ABSOLUTE)
  else()
    file(REAL_PATH ${prefix} prefix EXPAND_TILDE)
  endif()
  set(CMAKE_INSTALL_PREFIX ${prefix} CACHE PATH "top-level install path" FORCE)
endif()

set(LAPACK_ROOT ${CMAKE_INSTALL_PREFIX})
set(SCALAPACK_ROOT ${CMAKE_INSTALL_PREFIX})
set(MUMPS_ROOT ${CMAKE_INSTALL_PREFIX})
set(HDF5_ROOT ${CMAKE_INSTALL_PREFIX})
set(H5FORTRAN_ROOT ${CMAKE_INSTALL_PREFIX})
set(GLOW_ROOT ${CMAKE_INSTALL_PREFIX})

message(STATUS "Building and installing libraries under ${CMAKE_INSTALL_PREFIX}")
message(STATUS "please add system environment variables:")
foreach(p LAPACK_ROOT SCALAPACK_ROOT MUMPS_ROOT HDF5_ROOT H5FORTRAN_ROOT GLOW_ROOT)
  message(STATUS "   ${p}=${${p}}")
endforeach()

ExternalProject_Add(LAPACK
GIT_REPOSITORY ${lapack_git}
GIT_TAG ${lapack_tag}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${LAPACK_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
)

ExternalProject_Add(SCALAPACK
GIT_REPOSITORY ${scalapack_git}
GIT_TAG ${scalapack_tag}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${SCALAPACK_ROOT} -DLAPACK_ROOT:PATH=${LAPACK_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
DEPENDS LAPACK)

ExternalProject_Add(MUMPS
GIT_REPOSITORY ${mumps_git}
GIT_TAG ${mumps_tag}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${MUMPS_ROOT} -DLAPACK_ROOT:PATH=${LAPACK_ROOT} -DSCALAPACK_ROOT:PATH=${SCALAPACK_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
DEPENDS SCALAPACK)

include(${PROJECT_SOURCE_DIR}/../cmake/ext_libs/build_hdf5.cmake)

ExternalProject_Add(H5FORTRAN
GIT_REPOSITORY ${h5fortran_git}
GIT_TAG ${h5fortran_tag}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
CMAKE_ARGS -DHDF5_ROOT:PATH=${HDF5_ROOT} -DCMAKE_INSTALL_PREFIX:PATH=${H5FORTRAN_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false -Dhdf5_external=on
DEPENDS HDF5)

ExternalProject_Add(GLOW
GIT_REPOSITORY ${glow_git}
GIT_TAG ${glow_tag}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${GLOW_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
)
