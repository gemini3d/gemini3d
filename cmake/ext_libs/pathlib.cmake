include(ExternalProject)

find_package(FortranPathlib CONFIG QUIET)

if(FortranPathlib_FOUND)
  message(STATUS "Fortran Pathlib found: ${FortranPathlib_DIR}")
  return()
endif()

if(NOT FortranPathlib_ROOT)
  set(FortranPathlib_ROOT ${CMAKE_INSTALL_PREFIX})
endif()


if(BUILD_SHARED_LIBS)
  set(FortranPathlib_LIBRARIES ${FortranPathlib_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}pathlib${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(FortranPathlib_LIBRARIES ${FortranPathlib_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}pathlib${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(FortranPathlib_INCLUDE_DIRS ${FortranPathlib_ROOT}/include)

set(pathlib_cmake_args
--install-prefix=${FortranPathlib_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(PATHLIB
GIT_REPOSITORY ${pathlib_git}
GIT_TAG ${pathlib_tag}
CMAKE_ARGS ${pathlib_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${FortranPathlib_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${FortranPathlib_INCLUDE_DIRS})
# avoid generate race condition

add_library(FortranPathlib::pathlib INTERFACE IMPORTED)
target_link_libraries(FortranPathlib::pathlib INTERFACE ${FortranPathlib_LIBRARIES})
target_include_directories(FortranPathlib::pathlib INTERFACE ${FortranPathlib_INCLUDE_DIRS})

add_dependencies(FortranPathlib::pathlib PATHLIB)
