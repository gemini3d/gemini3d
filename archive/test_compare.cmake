function(octave_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})

add_test(NAME gemini:compare:${TESTNAME}:octave
COMMAND ${Octave_EXECUTABLE} --eval "compare_all('${_outdir}', '${_refdir}')")

set_tests_properties(gemini:compare:${TESTNAME}:octave PROPERTIES
TIMEOUT 30
REQUIRED_FILES ${_outdir}/inputs/config.nml
SKIP_RETURN_CODE 77
ENVIRONMENT OCTAVE_PATH=${PROJECT_SOURCE_DIR}/matlab
DISABLED ${octave_disabled})

endfunction(octave_compare)


function(matlab_compare TESTNAME)

set(_outdir ${PROJECT_BINARY_DIR}/test${TESTNAME})
set(_refdir ${PROJECT_SOURCE_DIR}/tests/data/test${TESTNAME})

add_test(NAME gemini:compare:${TESTNAME}:matlab
COMMAND ${Matlab_MAIN_PROGRAM} -batch "compare_all('${_outdir}', '${_refdir}')")

set_tests_properties(gemini:compare:${TESTNAME}:matlab PROPERTIES
TIMEOUT 60
REQUIRED_FILES ${_outdir}/inputs/config.nml
SKIP_RETURN_CODE 77
ENVIRONMENT MATLABPATH=${PROJECT_SOURCE_DIR}/matlab
DISABLED ${matlab_disabled})
# Matlab with a lot of toolboxes takes 15..30 seconds just to start, particularly on HPC with network file system.

endfunction(matlab_compare)