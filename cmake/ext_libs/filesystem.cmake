include(ExternalProject)

find_package(FortranFilesystem CONFIG QUIET)

if(FortranFilesystem_FOUND)
  message(STATUS "Fortran Filesystem found: ${FortranFilesystem_DIR}")
  return()
endif()

if(NOT FortranFilesystem_ROOT)
  set(FortranFilesystem_ROOT ${CMAKE_INSTALL_PREFIX})
endif()


if(BUILD_SHARED_LIBS)
  set(FortranFilesystem_LIBRARIES ${FortranFilesystem_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}filesystem${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(FortranFilesystem_LIBRARIES ${FortranFilesystem_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}filesystem${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(FortranFilesystem_INCLUDE_DIRS ${FortranFilesystem_ROOT}/include)

set(filesystem_cmake_args
--install-prefix=${FortranFilesystem_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(FortranFS
GIT_REPOSITORY ${filesystem_git}
GIT_TAG ${filesystem_tag}
CMAKE_ARGS ${filesystem_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${FortranFilesystem_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${FortranFilesystem_INCLUDE_DIRS})
# avoid generate race condition

add_library(FortranFilesystem::filesystem INTERFACE IMPORTED)
target_link_libraries(FortranFilesystem::filesystem INTERFACE ${FortranFilesystem_LIBRARIES})
target_include_directories(FortranFilesystem::filesystem INTERFACE ${FortranFilesystem_INCLUDE_DIRS})

add_dependencies(FortranFilesystem::filesystem FortranFS)
