# this enables CMake imported target HWM14::HWM14
include(ExternalProject)

# find_package(hwm14 CONFIG)

# if(hwm14_FOUND)
#   return()
# endif()

if(NOT HWM14_ROOT)
  set(HWM14_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if(BUILD_SHARED_LIBS)
  set(HWM14_LIBRARIES ${HWM14_ROOT}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}hwm14${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(HWM14_LIBRARIES ${HWM14_ROOT}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hwm14${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

set(hwm14_cmake_args
-DCMAKE_INSTALL_PREFIX:PATH=${HWM14_ROOT}
-DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
-DCMAKE_BUILD_TYPE=Release
-DBUILD_TESTING:BOOL=false
)

ExternalProject_Add(HWM14
GIT_REPOSITORY ${hwm14_git}
GIT_TAG ${hwm14_tag}
CMAKE_ARGS ${hwm14_cmake_args}
CMAKE_GENERATOR ${EXTPROJ_GENERATOR}
BUILD_BYPRODUCTS ${HWM14_LIBRARIES}
INACTIVITY_TIMEOUT 15
CONFIGURE_HANDLED_BY_BUILD ON
)

ExternalProject_Get_property(HWM14 SOURCE_DIR)


set(hwm14_dat_files
${SOURCE_DIR}/src/hwm14/hwm123114.bin
${SOURCE_DIR}/src/hwm14/dwm07b104i.dat
${SOURCE_DIR}/src/hwm14/gd2qd.dat
)
ExternalProject_Add_Step(HWM14 hwm_cp1 DEPENDEES update
COMMAND ${CMAKE_COMMAND} -E copy_if_different ${hwm14_dat_files} ${PROJECT_BINARY_DIR})
install(FILES ${hwm14_dat_files} DESTINATION bin)

add_library(HWM14::HWM14 INTERFACE IMPORTED)
target_link_libraries(HWM14::HWM14 INTERFACE "${HWM14_LIBRARIES}")

add_dependencies(HWM14::HWM14 HWM14)
