include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)

function(setup_gemini_test TESTNAME TIMEOUT)

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${TESTNAME})

if(python_ok AND hdf5)

  win32_hdf5_env()

  add_test(NAME gemini:hdf5:${TESTNAME}:dryrun
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_outdir} -dryrun)

  set_tests_properties(gemini:hdf5:${TESTNAME}:dryrun PROPERTIES
    TIMEOUT 60
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_SETUP hdf5:${TESTNAME}:dryrun)


  add_test(NAME gemini:hdf5:${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_outdir})

  # NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
  set_tests_properties(gemini:hdf5:${TESTNAME} PROPERTIES
    TIMEOUT ${TIMEOUT}
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_REQUIRED "hdf5:${TESTNAME}:dryrun;MPIMUMPS")

endif()

if(python_ok AND netcdf)
  add_test(NAME gemini:netcdf:${TESTNAME}:dryrun
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_outdir} -out_format nc -dryrun)

  set_tests_properties(gemini:netcdf:${TESTNAME}:dryrun PROPERTIES
    TIMEOUT 60
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_SETUP netcdf:${TESTNAME}:dryrun)

  add_test(NAME gemini:netcdf:${TESTNAME}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${TESTNAME} ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_outdir} -out_format nc)

  set_tests_properties(gemini:netcdf:${TESTNAME} PROPERTIES
    TIMEOUT ${TIMEOUT}
    SKIP_RETURN_CODE 77
    RUN_SERIAL true
    FIXTURES_REQUIRED "netcdf:${TESTNAME}:dryrun;MPIMUMPS")
endif()


compare_gemini_output(${TESTNAME})

endfunction(setup_gemini_test)
