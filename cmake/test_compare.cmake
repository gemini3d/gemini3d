function(python_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})


if(python_ok AND hdf5)
  add_test(NAME gemini:compare:hdf5:${TESTNAME}:python
  COMMAND gemini_compare ${_outdir} ${_refdir} -file_format h5)

  set_tests_properties(gemini:compare:hdf5:${TESTNAME}:python PROPERTIES
  TIMEOUT 30
  DEPENDS gemini:hdf5:${TESTNAME}
  REQUIRED_FILES ${_outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()


if(python_ok AND netcdf)
  add_test(NAME gemini:compare:netcdf:${TESTNAME}:python
  COMMAND gemini_compare ${_outdir} ${_refdir} -file_format nc)

  set_tests_properties(gemini:compare:netcdf:${TESTNAME}:python PROPERTIES
  TIMEOUT 30
  DEPENDS gemini:netcdf:${TESTNAME}
  REQUIRED_FILES ${_outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()

endfunction(python_compare)


function(compare_gemini_output TESTNAME)
python_compare(${TESTNAME})
endfunction(compare_gemini_output)
