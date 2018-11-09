cmake_minimum_required (VERSION 3.11)
project(interp12 Fortran)
enable_testing()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake/Modules/)

add_library(interp interpolation.f90)
set_target_properties(interp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/objects)

if(NOT TARGET const)
  add_library(const ${PROJECT_SOURCE_DIR}/../constants/phys_consts.f90)
endif()

target_link_libraries(interp PRIVATE const)
target_include_directories(interp PRIVATE ${CMAKE_BINARY_DIR}/numerical)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  include(test_interp.cmake)
endif()

