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


function(compare_gemini_output TESTNAME)
python_compare(${TESTNAME})
endfunction(compare_gemini_output)
