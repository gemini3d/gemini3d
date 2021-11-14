# for systems where you can't use "cmake -P scripts/install_prereq.cmake"
# this project builds the prereqs

cmake_minimum_required(VERSION 3.20...3.22)

project(BuildGeminiPrereq
LANGUAGES C Fortran
)

set(arith "s;d")
option(scotch "MUMPS scotch")
option(BUILD_SHARED_LIBS "Build shared libraries")
option(hdf5_parallel "HDF5 parallel")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  message(FATAL_ERROR "Specify where to install libraries by adding cmake option like:
   --install-prefix=/path/to/install")
endif()

# NOTE: lists must be passed via CMAKE_CACHE_ARGS using -Dvar:STRING=${var}
# the variable type must be included e.g. STRING
include(ExternalProject)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY EP_UPDATE_DISCONNECTED true)

include(${PROJECT_SOURCE_DIR}/../cmake/config/libraries.cmake)

set(HWLOC_ROOT ${CMAKE_INSTALL_PREFIX})
set(OPENMPI_ROOT ${CMAKE_INSTALL_PREFIX})
set(LAPACK_ROOT ${CMAKE_INSTALL_PREFIX})
set(SCALAPACK_ROOT ${CMAKE_INSTALL_PREFIX})
set(MUMPS_ROOT ${CMAKE_INSTALL_PREFIX})
set(ZLIB_ROOT ${CMAKE_INSTALL_PREFIX})
set(HDF5_ROOT ${CMAKE_INSTALL_PREFIX})
set(H5FORTRAN_ROOT ${CMAKE_INSTALL_PREFIX})
set(GLOW_ROOT ${CMAKE_INSTALL_PREFIX})

set(roots LAPACK_ROOT SCALAPACK_ROOT MUMPS_ROOT ZLIB_ROOT HDF5_ROOT H5FORTRAN_ROOT GLOW_ROOT)
if(NOT WIN32)
  list(APPEND roots HWLOC_ROOT OPENMPI_ROOT)
endif()
message(STATUS "Building and installing libraries under ${CMAKE_INSTALL_PREFIX}")
message(STATUS "please add system environment variables:")
foreach(p ${roots})
  message(STATUS "   ${p}=${${p}}")
endforeach()

# --- OpenMPI
if(NOT WIN32)

if(BUILD_SHARED_LIBS)
  set(hwloc_args --enable-shared --disable-static)
else()
  set(hwloc_args --disable-shared --enable-static)
endif()

find_package(LibXml2)
if(NOT LibXml2_FOUND)
  list(APPEND hwloc_args --disable-libxml2)
endif()

ExternalProject_Add(HWLOC
URL ${hwloc_url}
URL_HASH SHA256=${hwloc_sha256}
CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/HWLOC-prefix/src/HWLOC/configure --prefix=${HWLOC_ROOT} ${hwloc_args}
BUILD_COMMAND make -j
INSTALL_COMMAND make -j install
TEST_COMMAND ""
CONFIGURE_HANDLED_BY_BUILD ON
INACTIVITY_TIMEOUT 15
)

set(openmpi_flags
--prefix=${OPENMPI_ROOT}
--with-hwloc-libdir=${HWLOC_ROOT}/lib
)
if(BUILD_SHARED_LIBS)
  list(APPEND openmpi_flags --enable-shared --disable-static)
else()
  list(APPEND openmpi_flags --disable-shared --enable-static)
endif()
# --disable-shared avoids:
# lib_gcc/lib/libz.a(deflate_medium.c.o): relocation R_X86_64_32S against internal symbol `zng_length_code' can not be used when making a shared objec
# https://github.com/zlib-ng/zlib-ng/wiki/Common-build-problems#relocation-error-in-compress2

ExternalProject_Add(OPENMPI
URL ${openmpi_url}
URL_HASH SHA256=${openmpi_sha256}
CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/OPENMPI-prefix/src/OPENMPI/configure ${openmpi_flags}
BUILD_COMMAND make -j
INSTALL_COMMAND make -j install
TEST_COMMAND ""
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
DEPENDS HWLOC
)

endif()
# --- LAPACK
set(lapack_cmake_args
--install-prefix=${LAPACK_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(LAPACK
GIT_REPOSITORY ${lapack_git}
GIT_TAG ${lapack_tag}
CMAKE_ARGS ${lapack_cmake_args}
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
)

# --- Scalapack
set(scalapack_cmake_args
--install-prefix=${SCALAPACK_ROOT}
-DLAPACK_ROOT:PATH=${LAPACK_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dautobuild:BOOL=false
)

set(scalapack_deps LAPACK)
if(NOT WIN32)
  list(APPEND scalapack_deps OPENMPI)
endif()

ExternalProject_Add(SCALAPACK
GIT_REPOSITORY ${scalapack_git}
GIT_TAG ${scalapack_tag}
CMAKE_ARGS ${scalapack_cmake_args}
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
DEPENDS ${scalapack_deps}
)

# --- MUMPS
set(mumps_cmake_args
--install-prefix=${MUMPS_ROOT}
-DSCALAPACK_ROOT:PATH=${SCALAPACK_ROOT}
-DLAPACK_ROOT:PATH=${LAPACK_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dscotch:BOOL=${scotch}
-Dopenmp:BOOL=false
-Dparallel:BOOL=true
-Dautobuild:BOOL=false
)

ExternalProject_Add(MUMPS
GIT_REPOSITORY ${mumps_git}
GIT_TAG ${mumps_tag}
CMAKE_ARGS ${mumps_cmake_args}
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
DEPENDS SCALAPACK
)

# --- HDF5
set(zlib_cmake_args
-DZLIB_COMPAT:BOOL=on
-DZLIB_ENABLE_TESTS:BOOL=off
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
--install-prefix=${ZLIB_ROOT}
)

ExternalProject_Add(ZLIB
URL ${zlib_url}
URL_HASH SHA256=${zlib_sha256}
CMAKE_ARGS ${zlib_cmake_args}
CONFIGURE_HANDLED_BY_BUILD ON
INACTIVITY_TIMEOUT 15
)

set(hdf5_zlib
-DZLIB_ROOT:PATH=${ZLIB_ROOT}
-DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
-DZLIB_USE_EXTERNAL:BOOL=OFF
)

set(hdf5_cmake_args
${hdf5_zlib}
--install-prefix=${HDF5_ROOT}
-DCMAKE_MODULE_PATH:PATH=${CMAKE_MODULE_PATH}
-DHDF5_GENERATE_HEADERS:BOOL=false
-DHDF5_DISABLE_COMPILER_WARNINGS:BOOL=true
-DBUILD_STATIC_LIBS:BOOL=$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DHDF5_BUILD_FORTRAN:BOOL=true
-DHDF5_BUILD_CPP_LIB:BOOL=false
-DBUILD_TESTING:BOOL=false
-DHDF5_BUILD_EXAMPLES:BOOL=false
-DUSE_LIBAEC:bool=true
-DHDF5_BUILD_TOOLS:BOOL=$<NOT:$<BOOL:${hdf5_parallel}>>
-DHDF5_ENABLE_PARALLEL:BOOL=$<BOOL:${hdf5_parallel}>
-DMPI_ROOT:PATH=${MPI_ROOT}
)

ExternalProject_Add(HDF5
URL ${hdf5_url}
URL_HASH SHA256=${hdf5_sha256}
CMAKE_ARGS ${hdf5_cmake_args}
CONFIGURE_HANDLED_BY_BUILD ON
INACTIVITY_TIMEOUT 15
DEPENDS ZLIB
)

set(h5fortran_args
-DHDF5_ROOT:PATH=${HDF5_ROOT}
--install-prefix=${H5FORTRAN_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-Dautobuild:BOOL=off
)

set(h5fortran_deps HDF5)
if(hdf5_parallel AND NOT WIN32)
  list(APPEND h5fortran_deps OPENMPI)
endif()

ExternalProject_Add(H5FORTRAN
GIT_REPOSITORY ${h5fortran_git}
GIT_TAG ${h5fortran_tag}
CMAKE_ARGS ${h5fortran_args}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
DEPENDS ${h5fortran_deps}
)

set(glow_args
--install-prefix=${GLOW_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(GLOW
GIT_REPOSITORY ${glow_git}
GIT_TAG ${glow_tag}
CMAKE_ARGS ${glow_args}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
)
