# this enables CMake imported target HWM14::hwm14
include(ExternalProject)

if(NOT hwm14_external)
  find_package(hwm14 CONFIG)

  if(hwm14_FOUND)
    message(STATUS "HWM14 found: ${hwm14_DIR}")
    return()
  endif()
endif()

set(hwm14_external TRUE CACHE BOOL "Build HWM14")

if(NOT HWM14_ROOT)
  set(HWM14_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

set(HWM14_INCLUDE_DIRS ${HWM14_ROOT}/include)

if(BUILD_SHARED_LIBS)
  set(HWM14_LIBRARIES ${HWM14_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}hwm14${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(HWM14_LIBRARIES ${HWM14_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hwm14${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(hwm14_cmake_args
-DCMAKE_INSTALL_PREFIX=${HWM14_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(HWM14
GIT_REPOSITORY ${hwm14_git}
GIT_TAG ${hwm14_tag}
CMAKE_ARGS ${hwm14_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${HWM14_LIBRARIES}
TEST_COMMAND ""
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
)


set(hwm14_dat_files
${HWM14_ROOT}/share/data/HWM14/hwm123114.bin
${HWM14_ROOT}/share/data/HWM14/dwm07b104i.dat
${HWM14_ROOT}/share/data/HWM14/gd2qd.dat
)
ExternalProject_Add_Step(HWM14 hwm14_copy_data DEPENDEES install
COMMAND ${CMAKE_COMMAND} -E copy_if_different ${hwm14_dat_files} $<TARGET_FILE_DIR:gemini.bin>)
install(FILES ${hwm14_dat_files} TYPE BIN)


add_library(HWM14::hwm14 INTERFACE IMPORTED)
target_link_libraries(HWM14::hwm14 INTERFACE "${HWM14_LIBRARIES}")
target_include_directories(HWM14::hwm14 INTERFACE ${HWM14_INCLUDE_DIRS})

add_dependencies(HWM14::hwm14 HWM14)
