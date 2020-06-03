function(python_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})


if(python_ok AND hdf5)
  add_test(NAME gemini:compare:hdf5:${TESTNAME}:python
  COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/compare_all.py ${_outdir} ${_refdir} -file_format h5)

  set_tests_properties(gemini:compare:hdf5:${TESTNAME}:python PROPERTIES
  TIMEOUT 30
  DEPENDS gemini:hdf5:${TESTNAME}
  REQUIRED_FILES ${_outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()


if(python_ok AND netcdf)
  add_test(NAME gemini:compare:netcdf:${TESTNAME}:python
  COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/compare_all.py ${_outdir} ${_refdir} -file_format nc)

  set_tests_properties(gemini:compare:netcdf:${TESTNAME}:python PROPERTIES
  TIMEOUT 30
  DEPENDS gemini:netcdf:${TESTNAME}
  REQUIRED_FILES ${_outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()

endfunction(python_compare)


function(h5diff_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})

set(_reltol 0.15)
# 0.15 <==> 15% relative tolerance
# this large relative tolerance is because h5diff can only do rel or abs, while Python can do both
# and so h5diff has to be looser for miniscule numerical noise.

if(HDF5_DIFF_EXECUTABLE)
  add_test(NAME gemini:compare:${TESTNAME}:h5diff
  COMMAND ${HDF5_DIFF_EXECUTABLE} --relative=${_reltol} ${_outdir}/20130220_18300.000000.h5 ${_refdir}/20130220_18300.000000.h5)

  set_tests_properties(gemini:compare:${TESTNAME}:h5diff PROPERTIES
  DEPENDS gemini:hdf5:${TESTNAME}
  TIMEOUT 15)
endif()

endfunction(h5diff_compare)


function(compare_gemini_output TESTNAME)
h5diff_compare(${TESTNAME})
python_compare(${TESTNAME})
endfunction(compare_gemini_output)
