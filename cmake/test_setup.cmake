include(${CMAKE_CURRENT_LIST_DIR}/test_compare.cmake)

function(setup_gemini_test testname TIMEOUT)

# --- setup test
set(_outdir ${CMAKE_CURRENT_BINARY_DIR}/test${testname})

if(python_ok)

if(mpi)
  set(_cmd ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${testname} -mpiexec ${MPIEXEC_EXECUTABLE} $<TARGET_FILE:gemini.bin> ${_outdir})
else()
  set(_cmd ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_test.py ${testname} $<TARGET_FILE:gemini.bin> ${_outdir})
endif(mpi)

if(hdf5)

add_test(NAME gemini:hdf5:${testname}:dryrun
  COMMAND ${_cmd} -dryrun
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  # NOTE: Working_Diretory is NECESSARY for Windows + Intel + HDF5

set_tests_properties(gemini:hdf5:${testname}:dryrun PROPERTIES
  TIMEOUT 60
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  FIXTURES_REQUIRED mumps_fixture
  FIXTURES_SETUP hdf5:${testname}:dryrun)


add_test(NAME gemini:hdf5:${testname}
  COMMAND ${_cmd}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  # NOTE: Working_Directory is NECESSARY for Windows + Intel + HDF5

# NOTE: don't use REQUIRED_FILES because it won't let file download if not present.
set_tests_properties(gemini:hdf5:${testname} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  FIXTURES_REQUIRED hdf5:${testname}:dryrun
  FIXTURES_SETUP hdf5:${testname})

endif(hdf5)

if(netcdf)
add_test(NAME gemini:netcdf:${testname}:dryrun
  COMMAND ${_cmd} -out_format nc -dryrun)

set_tests_properties(gemini:netcdf:${testname}:dryrun PROPERTIES
  TIMEOUT 60
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  FIXTURES_REQUIRED mumps_fixture
  FIXTURES_SETUP netcdf:${testname}:dryrun)

add_test(NAME gemini:netcdf:${testname}
  COMMAND ${_cmd} -out_format nc)

set_tests_properties(gemini:netcdf:${testname} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  RUN_SERIAL true
  FIXTURES_REQUIRED netcdf:${testname}:dryrun
  FIXTURES_SETUP netcdf:${testname})
endif(netcdf)

endif(python_ok)

compare_gemini_output(${testname})

endfunction(setup_gemini_test)
