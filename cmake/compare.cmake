include(${CMAKE_CURRENT_LIST_DIR}/python.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/octave.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/matlab.cmake)

function(python_compare TESTNAME OUTDIR REFDIR)

add_test(NAME gemini:compare:${TESTNAME}:python
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_all.py ${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR} ${REFDIR})

set_tests_properties(gemini:compare:${TESTNAME}:python PROPERTIES
  TIMEOUT 30
  REQUIRED_FILES ${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/inputs/config.nml
  SKIP_RETURN_CODE 77)

endfunction(python_compare)


function(octave_compare TESTNAME OUTDIR REFDIR)

add_test(NAME gemini:compare:${TESTNAME}:octave
  COMMAND ${Octave_EXECUTABLE} --eval "compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${REFDIR}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(gemini:compare:${TESTNAME}:octave PROPERTIES
  TIMEOUT 30
  REQUIRED_FILES ${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/inputs/config.nml
  SKIP_RETURN_CODE 77)

endfunction(octave_compare)


function(matlab_compare TESTNAME OUTDIR REFDIR)

add_test(NAME gemini:compare:${TESTNAME}:matlab
  COMMAND ${Matlab_MAIN_PROGRAM} -batch compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${REFDIR}')
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(gemini:compare:${TESTNAME}:matlab PROPERTIES
  TIMEOUT 60
  REQUIRED_FILES ${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/inputs/config.nml
  SKIP_RETURN_CODE 77)
# Matlab with a lot of toolboxes takes 15..30 seconds just to start, particularly on HPC with network file system.

endfunction(matlab_compare)


function(compare_gemini_output TESTNAME OUTDIR REFDIR)
# This sets up the Compare* tests

if(PythonOK)
  python_compare(${TESTNAME} ${OUTDIR} ${REFDIR})
endif(PythonOK)

if(OctaveOK)
  octave_compare(${TESTNAME} ${OUTDIR} ${REFDIR})
endif(OctaveOK)

if(MatlabOK)
  matlab_compare(${TESTNAME} ${TESTDIR} ${REFDIR})
endif(MatlabOK)

endfunction(compare_gemini_output)