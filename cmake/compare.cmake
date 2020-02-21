include(${CMAKE_CURRENT_LIST_DIR}/python.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/octave.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/matlab.cmake)

function(python_compare TESTNAME OUTDIR REFDIR)

add_test(NAME ${TESTNAME}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/compare_all.py ${CMAKE_BINARY_DIR}/${OUTDIR} ${REFDIR})

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 30
  SKIP_RETURN_CODE 77)

endfunction(python_compare)


function(octave_compare TESTNAME OUTDIR REFDIR)

add_test(NAME ${TESTNAME}
  COMMAND ${Octave_EXECUTABLE} --eval "compare_all('${CMAKE_BINARY_DIR}/${OUTDIR}','${REFDIR}')"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 30
  SKIP_RETURN_CODE 77)

endfunction(octave_compare)


function(matlab_compare TESTNAME OUTDIR REFDIR)

add_test(NAME ${TESTNAME}
  COMMAND ${Matlab_MAIN_PROGRAM} -batch compare_all('${CMAKE_BINARY_DIR}/${OUTDIR}','${REFDIR}')
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 60
  SKIP_RETURN_CODE 77)
# Matlab with a lot of toolboxes takes 15..30 seconds just to start, particularly on HPC with network file system.

endfunction(matlab_compare)


function(compare_gemini_output TESTNAME OUTDIR REFDIR)
# This sets up the Compare* tests

if(PythonOK)
  python_compare(${TESTNAME}:python ${OUTDIR} ${REFDIR})
endif(PythonOK)

if(OctaveOK)
  octave_compare(${TESTNAME}:octave ${OUTDIR} ${REFDIR})
endif(OctaveOK)

if(MatlabOK)
  matlab_compare(${TESTNAME}:matlab ${TESTDIR} ${REFDIR})
endif(MatlabOK)

endfunction(compare_gemini_output)