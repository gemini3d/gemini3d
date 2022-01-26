include(ExternalProject)

if(NOT ffilesystem_external)
  find_package(ffilesystem CONFIG)

  if(ffilesystem_FOUND)
    message(STATUS "Fortran Filesystem found: ${ffilesystem_DIR}")
    return()
  endif()
endif()

set(ffilesystem_external true CACHE BOOL "build Fortran-filesystem")

if(NOT ffilesystem_ROOT)
  set(ffilesystem_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

set(ffs_implib)
if(BUILD_SHARED_LIBS)
  if(WIN32)
    set(ffs_implib ${ffilesystem_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}filesystem${CMAKE_SHARED_LIBRARY_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(ffilesystem_LIBRARIES ${ffilesystem_ROOT}/bin/${CMAKE_SHARED_LIBRARY_PREFIX}filesystem${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(ffilesystem_LIBRARIES ${ffilesystem_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}filesystem${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()
else()
  set(ffilesystem_LIBRARIES ${ffilesystem_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}filesystem${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(ffilesystem_INCLUDE_DIRS ${ffilesystem_ROOT}/include)

set(ffilesystem_cmake_args
-DCMAKE_INSTALL_PREFIX=${ffilesystem_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
)

ExternalProject_Add(FFILESYSTEM
GIT_REPOSITORY ${ffilesystem_git}
GIT_TAG ${ffilesystem_tag}
CMAKE_ARGS ${ffilesystem_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${ffilesystem_LIBRARIES} ${ffs_implib}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD true
)

file(MAKE_DIRECTORY ${ffilesystem_INCLUDE_DIRS})
# avoid generate race condition

if(BUILD_SHARED_LIBS)
  add_library(ffilesystem::filesystem SHARED IMPORTED)
  if(WIN32)
    set_target_properties(ffilesystem::filesystem PROPERTIES IMPORTED_IMPLIB ${ffs_implib})
  endif()
else()
  add_library(ffilesystem::filesystem STATIC IMPORTED)
endif()

set_target_properties(ffilesystem::filesystem PROPERTIES
IMPORTED_LOCATION ${ffilesystem_LIBRARIES}
INTERFACE_INCLUDE_DIRECTORIES ${ffilesystem_INCLUDE_DIRS}
LINKER_LANGUAGE CXX
)
target_link_libraries(ffilesystem::filesystem INTERFACE ${lib_filesystem})

# target_link_libraries(ffilesystem::filesystem INTERFACE stdc++)  # did not help
# instead, set linker_langauge CXX for the specific targets linking ffilesystem::filesystem

add_dependencies(ffilesystem::filesystem FFILESYSTEM)
