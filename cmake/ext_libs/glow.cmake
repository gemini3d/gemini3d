include(ExternalProject)

if(NOT glow_external)
  find_package(glow CONFIG)

  if(glow_FOUND)
    message(STATUS "GLOW found: ${glow_DIR}")
    return()
  endif()
endif()

set(glow_external true CACHE BOOL "build GLOW")

if(NOT GLOW_ROOT)
  set(GLOW_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(GLOW_LIBRARIES
  ${GLOW_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}glow${CMAKE_SHARED_LIBRARY_SUFFIX}
  ${GLOW_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}msis00${CMAKE_SHARED_LIBRARY_SUFFIX}
  )
else()
  set(GLOW_LIBRARIES
  ${GLOW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}glow${CMAKE_STATIC_LIBRARY_SUFFIX}
  ${GLOW_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}msis00${CMAKE_STATIC_LIBRARY_SUFFIX}
  )
endif()

set(GLOW_INCLUDE_DIRS ${GLOW_ROOT}/include)

set(glow_cmake_args
-DCMAKE_INSTALL_PREFIX=${GLOW_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(GLOW
GIT_REPOSITORY ${glow_git}
GIT_TAG ${glow_tag}
CMAKE_ARGS ${glow_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${GLOW_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${GLOW_INCLUDE_DIRS})
# avoid generate race condition

add_library(glow::glow INTERFACE IMPORTED)
target_link_libraries(glow::glow INTERFACE ${GLOW_LIBRARIES})
target_include_directories(glow::glow INTERFACE ${GLOW_INCLUDE_DIRS})
target_compile_definitions(glow::glow INTERFACE DATADIR="${GLOW_ROOT}/share/data/glow/")
# DATADIR comes from glow project BUILD_INTERFACE COMPILE_DEFINITIONS
# this does not seem to work with get_target_property()
# so preprocess instead of configure_file()

add_dependencies(glow::glow GLOW)
