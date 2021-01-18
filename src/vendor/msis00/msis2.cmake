include(FetchContent)

if(NOT DEFINED msis2_url OR NOT DEFINED msis2_sha1)
  if(CMAKE_VERSION VERSION_LESS 3.19)
    message(FATAL_ERROR "CMake >= 3.19 required to use MSIS 2.0")
  endif()

  set(_json_file ${PROJECT_SOURCE_DIR}/cmake/libraries.json)
  if(NOT EXISTS ${_json_file})
    message(FATAL_ERROR "must define msis2_url and msis2_sha1")
  endif()
  file(READ ${_json_file} _libj)

  string(JSON msis2_url GET ${_libj} "msis2" "url")
  string(JSON msis2_sha1 GET ${_libj} "msis2" "sha1")
endif()

FetchContent_Declare(msis2proj
URL ${msis2_url}
URL_HASH SHA1=${msis2_sha1}
UPDATE_DISCONNECTED true
)

FetchContent_MakeAvailable(msis2proj)

set(_s ${msis2proj_SOURCE_DIR})  # convenience

add_library(msis2 ${_s}/alt2gph.F90 ${_s}/msis_constants.F90 ${_s}/msis_init.F90 ${_s}/msis_gfn.F90 ${_s}/msis_tfn.F90 ${_s}/msis_dfn.F90 ${_s}/msis_calc.F90 ${_s}/msis_gtd8d.F90)
set_target_properties(msis2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(msis2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # msis_calc:bspline has argument mismatch on nodes variable
  target_compile_options(msis2 PRIVATE -std=legacy)
endif()

# MSIS 2.0 needs this parm file.
# From your Fortran code, refer to this file by
# `call msisinit(parmpath=)` perhaps via CMake configure_file()
file(COPY ${msis2proj_SOURCE_DIR}/msis20.parm DESTINATION ${PROJECT_BINARY_DIR})

if(BUILD_TESTING)
  add_executable(msis2test ${msis2proj_SOURCE_DIR}/msis2.0_test.F90)
  target_link_libraries(msis2test PRIVATE msis2)

  add_test(NAME MSIS2
    COMMAND $<TARGET_FILE:msis2test>
    WORKING_DIRECTORY ${msis2proj_SOURCE_DIR})
endif()
