include(FetchContent)

FetchContent_Declare(MSIS2
URL ${msis2_zip}
URL_HASH SHA256=${msis2_sha256}
INACTIVITY_TIMEOUT 15
)

FetchContent_MakeAvailable(MSIS2)

set(msis2_SOURCE_DIR ${msis2_SOURCE_DIR} PARENT_SCOPE)
set(_s ${msis2_SOURCE_DIR})
# convenience

# patching API MSIS
add_custom_command(
OUTPUT ${msis2_BINARY_DIR}/msis_calc.F90
COMMAND ${CMAKE_COMMAND} -Din_file:FILEPATH=${msis2_SOURCE_DIR}/msis_calc.F90 -Dpatch_file:FILEPATH=${PROJECT_SOURCE_DIR}/src/vendor/nrl_msis/msis_api.patch -Dout_file:FILEPATH=${msis2_BINARY_DIR}/msis_calc.F90 -P ${PROJECT_SOURCE_DIR}/cmake/PatchFile.cmake
DEPENDS ${msis2_SOURCE_DIR}/msis_calc.F90
)

add_library(msis2mod
${_s}/alt2gph.F90
${_s}/msis_constants.F90
${_s}/msis_init.F90
${_s}/msis_gfn.F90
${_s}/msis_tfn.F90
${_s}/msis_dfn.F90
${_s}/msis_gtd8d.F90
${msis2_BINARY_DIR}/msis_calc.F90
)

# MSIS 2.0 needs this parm file.
# NOTE: there is also add_custom_command() to copy this to each config's binary dir in src/cmakelists.txt
install(FILES ${msis2_SOURCE_DIR}/msis20.parm TYPE BIN)

if(${PROJECT}_BUILD_TESTING)
  add_executable(msis2test ${msis2_SOURCE_DIR}/msis2.0_test.F90)
  target_link_libraries(msis2test PRIVATE msis2)

  add_test(NAME MSIS2
  COMMAND $<TARGET_FILE:msis2test>
  WORKING_DIRECTORY ${msis2_SOURCE_DIR}
  )
endif()
