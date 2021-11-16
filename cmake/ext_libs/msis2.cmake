include(FetchContent)

FetchContent_Declare(MSIS2
URL ${msis2_zip}
URL_HASH SHA256=${msis2_sha256}
INACTIVITY_TIMEOUT 15
)

FetchContent_MakeAvailable(MSIS2)

set(_s ${msis2_SOURCE_DIR})  # convenience

add_library(msis2mod
${_s}/alt2gph.F90
${_s}/msis_constants.F90
${_s}/msis_init.F90
${_s}/msis_gfn.F90
${_s}/msis_tfn.F90
${_s}/msis_dfn.F90
${_s}/msis_calc.F90
${_s}/msis_gtd8d.F90
)

# MSIS 2.0 needs this parm file.
file(COPY ${msis2_SOURCE_DIR}/msis20.parm DESTINATION ${PROJECT_BINARY_DIR})
install(FILES ${msis2_SOURCE_DIR}/msis20.parm TYPE BIN)

if(${PROJECT}_BUILD_TESTING)
  add_executable(msis2test ${msis2_SOURCE_DIR}/msis2.0_test.F90)
  target_link_libraries(msis2test PRIVATE msis2)

  add_test(NAME MSIS2
  COMMAND $<TARGET_FILE:msis2test>
  WORKING_DIRECTORY ${msis2_SOURCE_DIR}
  )
endif()

# patching API MSIS

if(msis_patched)
  return()
endif()

set(msis_orig ${msis2_SOURCE_DIR}/msis_calc.F90)
set(msis_patch ${PROJECT_SOURCE_DIR}/src/vendor/nrl_msis/msis_api.patch)

find_package(PATCH)

patch_file(${msis_orig} ${msis_patch})

set(msis_patched true CACHE BOOL "MSIS is patched")
