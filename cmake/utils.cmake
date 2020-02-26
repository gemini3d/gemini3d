include(${CMAKE_CURRENT_LIST_DIR}/compare.cmake)

function(setup_gemini_test TESTNAME EXE TESTDIR REFDIR TIMEOUT)

if(NOT PythonOK)
  message(VERBOSE " SKIP ${TESTNAME}: Python not found.")
  return()
endif()

set(_config_file ${CMAKE_CURRENT_SOURCE_DIR}/initialize/${TESTDIR}/config.nml)

add_test(NAME gemini:${TESTNAME}
  COMMAND ${Python3_EXECUTABLE} script_utils/meson_run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:${EXE}> ${_config_file} ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(gemini:${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  REQUIRED_FILES ${_config_file}
  FIXTURES_REQUIRED MPIMUMPS
  FIXTURES_SETUP GeminiRun_${TESTNAME}
  RUN_SERIAL true
)
if(WIN32 AND HDF5_ROOT)  # for Windows ifort dll
  set_tests_properties(gemini:${TESTNAME} PROPERTIES ENVIRONMENT "PATH=${HDF5_ROOT}/bin;$ENV{PATH}")
endif()

endfunction(setup_gemini_test)
