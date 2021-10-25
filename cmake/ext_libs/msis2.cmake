include(FetchContent)

FetchContent_Declare(MSIS2
URL ${msis2_zip}
URL_HASH SHA256=${msis2_sha256}
INACTIVITY_TIMEOUT 15
)

FetchContent_MakeAvailable(MSIS2)

set(_s ${msis2_SOURCE_DIR})  # convenience

add_library(msis2
${_s}/alt2gph.F90
${_s}/msis_constants.F90
${_s}/msis_init.F90
${_s}/msis_gfn.F90
${_s}/msis_tfn.F90
${_s}/msis_dfn.F90
${_s}/msis_calc.F90
${_s}/msis_gtd8d.F90)

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

# patching API MSIS

if(msis_patched)
  return()
endif()

set(msis_orig ${msis2_SOURCE_DIR}/msis_calc.F90)
set(msis_patch ${PROJECT_SOURCE_DIR}/src/vendor/nrl_msis/msis_api.patch)
# API patch
if(WIN32)
  find_program(WSL NAMES wsl REQUIRED)

  execute_process(COMMAND ${WSL} wslpath ${msis_orig}
    TIMEOUT 5
    OUTPUT_VARIABLE msis_orig_path
    COMMAND_ERROR_IS_FATAL ANY
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(COMMAND ${WSL} wslpath ${msis_patch}
    TIMEOUT 5
    OUTPUT_VARIABLE msis_patch_path
    COMMAND_ERROR_IS_FATAL ANY
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(COMMAND ${WSL} patch ${msis_orig_path} ${msis_patch_path}
    TIMEOUT 10
    COMMAND_ERROR_IS_FATAL ANY
  )
else()
  find_program(PATCH NAMES patch REQUIRED)
  execute_process(COMMAND ${PATCH} ${msis_orig} ${msis_patch}
    TIMEOUT 10
    COMMAND_ERROR_IS_FATAL ANY
  )
endif()

set(msis_patched true CACHE BOOL "MSIS is patched")
