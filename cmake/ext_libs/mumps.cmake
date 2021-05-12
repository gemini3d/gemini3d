# Intel MKL-compiled MUMPS requires at the linker for the main executable:
# mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core
#
# easily obtain MUMPS without compiling:
# CentOS 6/7 EPEL: yum install mumps-devel
# Ubuntu / Debian: apt install libmumps-dev

include(ExternalProject)

# --- prereqs
include(${CMAKE_CURRENT_LIST_DIR}/lapack.cmake)

if(mpi)
  include(${CMAKE_CURRENT_LIST_DIR}/scalapack.cmake)

  find_package(HWLOC)  # need here for cmake/summary.cmake scope
endif()
# --- MUMPS

if(NOT mumps_external AND (MUMPS_ROOT OR (DEFINED ENV{MUMPS_ROOT}) OR (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)))
  set(_comp ${arith})
  if(NOT mpi)
    list(APPEND _comp mpiseq)
  endif()
  if(scotch)
    list(APPEND _comp Scotch)
  endif()
  if(metis)
    list(APPEND _comp METIS)
  endif()
  if(OpenMP_FOUND)
    list(APPEND _comp OpenMP)
  endif()

  if(autobuild)
    find_package(MUMPS COMPONENTS ${_comp})
  else()
    find_package(MUMPS COMPONENTS ${_comp} REQUIRED)
  endif()
endif()

if(MUMPS_FOUND OR TARGET MUMPS::MUMPS)
  return()
endif()

set(mumps_external true CACHE BOOL "build Mumps")

if(SCALAPACK_FOUND OR NOT scalapack_external)
  add_custom_target(SCALAPACK)
endif()

if(NOT MUMPS_ROOT)
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(MUMPS_ROOT ${PROJECT_BINARY_DIR})
  else()
    set(MUMPS_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()
endif()

set(MUMPS_INCLUDE_DIRS ${MUMPS_ROOT}/include)
set(MUMPS_LIBRARIES)

foreach(a ${arith})
  list(APPEND MUMPS_LIBRARIES ${MUMPS_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${a}mumps${CMAKE_STATIC_LIBRARY_SUFFIX})
endforeach()

list(APPEND MUMPS_LIBRARIES
${MUMPS_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mumps_common${CMAKE_STATIC_LIBRARY_SUFFIX}
${MUMPS_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}pord${CMAKE_STATIC_LIBRARY_SUFFIX})

if(NOT MPI_FOUND)
  set(MUMPS_MPISEQ_LIBRARIES ${MUMPS_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${a}mpiseq${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

ExternalProject_Add(MUMPS
GIT_REPOSITORY ${mumps_git}
GIT_TAG ${mumps_tag}
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${MUMPS_ROOT} -DSCALAPACK_ROOT:PATH=${SCALAPACK_ROOT} -DLAPACK_ROOT:PATH=${LAPACK_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false -Dmetis:BOOL=${metis} -Dscotch:BOOL=${scotch} -Dopenmp:BOOL=false -Dparallel:BOOL=${mpi}
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
BUILD_BYPRODUCTS ${MUMPS_LIBRARIES} ${MUMPS_MPISEQ_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
DEPENDS SCALAPACK
)

file(MAKE_DIRECTORY ${MUMPS_INCLUDE_DIRS})

add_library(MUMPS::MUMPS INTERFACE IMPORTED)
target_link_libraries(MUMPS::MUMPS INTERFACE "${MUMPS_LIBRARIES}")
target_include_directories(MUMPS::MUMPS INTERFACE ${MUMPS_INCLUDE_DIRS})

# race condition for linking without this
add_dependencies(MUMPS::MUMPS MUMPS)

if(NOT MPI_FOUND)
  add_library(MUMPS::MPISEQ INTERFACE IMPORTED)
  target_link_libraries(MUMPS::MPISEQ INTERFACE "${MUMPS_MPISEQ_LIBRARIES}")
  target_include_directories(MUMPS::MPISEQ INTERFACE ${MUMPS_INCLUDE_DIRS})

  # race condition for linking without this
  add_dependencies(MUMPS::MPISEQ MUMPS)
endif()
