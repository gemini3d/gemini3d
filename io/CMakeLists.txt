cmake_minimum_required(VERSION 3.1)
project(filesystemutils Fortran)
enable_testing()

get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent)
  if(NOT realbits)
    set(realbits 64)
  endif()

  include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/compilers.cmake)

  add_library(const ${CMAKE_CURRENT_SOURCE_DIR}/../numerical/constants/phys_consts.F90)
  target_compile_definitions(const PRIVATE REALBITS=${realbits})
  set_target_properties(const PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

  add_library(timeutils ${CMAKE_CURRENT_SOURCE_DIR}/../temporal/timeutils.f90)
  target_link_libraries(timeutils PRIVATE const)
  set_target_properties(timeutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

add_library(fsutils mkdir.F90 expanduser.f90 formats.f90)
target_compile_options(fsutils PRIVATE ${FFLAGS})
set_target_properties(fsutils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(logging logging.f90)
target_compile_options(logging PRIVATE ${FFLAGS})
set_target_properties(logging PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_formats test_formats.f90)
target_compile_options(test_formats PRIVATE ${FFLAGS})
target_link_libraries(test_formats PRIVATE fsutils timeutils)
set_target_properties(test_formats PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME DateFormats COMMAND test_formats)
set_tests_properties(DateFormats PROPERTIES TIMEOUT 5.0
  FIXTURES_SETUP IOfmt)
