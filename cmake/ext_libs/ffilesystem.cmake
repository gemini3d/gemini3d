include(ExternalProject)

find_package(ffilesystem CONFIG QUIET)

if(ffilesystem_FOUND)
  message(STATUS "Fortran Filesystem found: ${ffilesystem_DIR}")
  return()
endif()

if(NOT ffilesystem_ROOT)
  set(ffilesystem_ROOT ${CMAKE_INSTALL_PREFIX})
endif()


if(BUILD_SHARED_LIBS)
  set(ffilesystem_LIBRARIES ${ffilesystem_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}filesystem${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(ffilesystem_LIBRARIES ${ffilesystem_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}filesystem${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(ffilesystem_INCLUDE_DIRS ${ffilesystem_ROOT}/include)

set(ffilesystem_cmake_args
--install-prefix=${ffilesystem_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(FFILESYSTEM
GIT_REPOSITORY ${ffilesystem_git}
GIT_TAG ${ffilesystem_tag}
CMAKE_ARGS ${ffilesystem_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${ffilesystem_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${ffilesystem_INCLUDE_DIRS})
# avoid generate race condition

add_library(ffilesystem::filesystem INTERFACE IMPORTED)
target_link_libraries(ffilesystem::filesystem INTERFACE ${ffilesystem_LIBRARIES})
target_include_directories(ffilesystem::filesystem INTERFACE ${ffilesystem_INCLUDE_DIRS})

add_dependencies(ffilesystem::filesystem FFILESYSTEM)
