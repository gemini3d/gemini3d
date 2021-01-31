include(FetchContent)

if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "MSIS 2.0 requires CMake >= 3.19")
endif()

FetchContent_Declare(MSIS2
  URL ${msis2_zip}
  URL_HASH SHA1=${msis2_sha1}
  TLS_VERIFY ON)
FetchContent_MakeAvailable(MSIS2)

set(_s ${msis2_SOURCE_DIR})  # convenience

add_library(msis2 ${_s}/alt2gph.F90 ${_s}/msis_constants.F90 ${_s}/msis_init.F90 ${_s}/msis_gfn.F90 ${_s}/msis_tfn.F90 ${_s}/msis_dfn.F90 ${_s}/msis_calc.F90 ${_s}/msis_gtd8d.F90)
set_target_properties(msis2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(msis2 INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)

# MSIS 2.0 needs this parm file.
file(COPY ${msis2_SOURCE_DIR}/msis20.parm DESTINATION ${PROJECT_BINARY_DIR})
install(FILES ${msis2_SOURCE_DIR}/msis20.parm DESTINATION bin)

if(${PROJECT}_BUILD_TESTING)
  add_executable(msis2test ${msis2_SOURCE_DIR}/msis2.0_test.F90)
  target_link_libraries(msis2test PRIVATE msis2)

  add_test(NAME MSIS2
    COMMAND $<TARGET_FILE:msis2test>
    WORKING_DIRECTORY ${msis2_SOURCE_DIR})
endif()
