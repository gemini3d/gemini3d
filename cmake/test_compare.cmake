include(${CMAKE_CURRENT_LIST_DIR}/python.cmake)


function(python_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})

add_test(NAME gemini:compare:${TESTNAME}:python
COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/compare_all.py ${_outdir} ${_refdir})

set_tests_properties(gemini:compare:${TESTNAME}:python PROPERTIES
TIMEOUT 30
REQUIRED_FILES ${_outdir}/inputs/config.nml
SKIP_RETURN_CODE 77
DISABLED ${python_disabled})

endfunction(python_compare)


function(h5diff_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})

set(_reltol 0.01)  # 0.01 <==> 1% relative tolerance

add_test(NAME gemini:compare:${TESTNAME}:h5diff
COMMAND ${HDF5_DIFF_EXECUTABLE} --relative=${_reltol} ${_outdir}/20130220_18300.000000.h5 ${_refdir}/20130220_18300.000000.h5)

set(h5diff_disabled false)
if(NOT HDF5_DIFF_EXECUTABLE)
  set(h5diff_disabled true)
endif()

set_tests_properties(gemini:compare:${TESTNAME}:h5diff PROPERTIES
TIMEOUT 15
REQUIRED_FILES "${_outdir}/20130220_18300.000000.h5;${_refdir}/20130220_18300.000000.h5"
DISABLED ${h5diff_disabled})

endfunction(h5diff_compare)


function(compare_gemini_output TESTNAME)
h5diff_compare(${TESTNAME})
python_compare(${TESTNAME})
endfunction(compare_gemini_output)
