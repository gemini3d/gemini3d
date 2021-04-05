# Finds Lapack, tests, and if not found or broken, autobuild Lapack
include(ExternalProject)

if(NOT lapack_external)
  if(autobuild)
    find_package(LAPACK)
  else()
    find_package(LAPACK REQUIRED)
  endif()
endif()

if(LAPACK_FOUND OR TARGET LAPACK::LAPACK)
  return()
endif()


set(lapack_external true CACHE BOOL "build Lapack")

if(NOT DEFINED LAPACK_ROOT)
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(LAPACK_ROOT ${PROJECT_BINARY_DIR}/lapack)
  else()
    set(LAPACK_ROOT ${CMAKE_INSTALL_PREFIX})
  endif()
endif()

set(LAPACK_LIBRARIES
${LAPACK_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}
${LAPACK_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX})


ExternalProject_Add(LAPACK
GIT_REPOSITORY ${lapack_git}
GIT_TAG ${lapack_tag}
CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${LAPACK_ROOT} -DBUILD_SHARED_LIBS:BOOL=false -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING:BOOL=false
CMAKE_CACHE_ARGS -Darith:STRING=${arith}
BUILD_BYPRODUCTS ${LAPACK_LIBRARIES}
INACTIVITY_TIMEOUT 30
CONFIGURE_HANDLED_BY_BUILD ON
)

add_library(LAPACK::LAPACK INTERFACE IMPORTED)
target_link_libraries(LAPACK::LAPACK INTERFACE "${LAPACK_LIBRARIES}")

# race condition for linking without this
add_dependencies(LAPACK::LAPACK LAPACK)
