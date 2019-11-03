if(usehdf5)
  set(OUTEXT ".h5")
else()
  set(OUTEXT ".dat")
endif()

set(_zenodo_ext ".dat")

function(python_compare TESTNAME OUTDIR REFDIR REQSTEM)

add_test(NAME ${TESTNAME}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/compare_all.py
    ${CMAKE_BINARY_DIR}/${OUTDIR} ${CMAKE_SOURCE_DIR}/${REFDIR})

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 30
  REQUIRED_FILES "${CMAKE_BINARY_DIR}/${OUTDIR}/${REQSTEM}${OUTEXT};${CMAKE_SOURCE_DIR}/${REFDIR}/${REQSTEM}${_zenodo_ext}")

endfunction(python_compare)


function(octave_compare TESTNAME OUTDIR REFDIR REQSTEM)

add_test(NAME ${TESTNAME}
  COMMAND ${Octave_EXECUTABLE} --eval "compare_all('${CMAKE_BINARY_DIR}/${OUTDIR}','${CMAKE_SOURCE_DIR}/${REFDIR}')"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 30
  REQUIRED_FILES "${CMAKE_BINARY_DIR}/${OUTDIR}/${REQSTEM}${OUTEXT};${CMAKE_SOURCE_DIR}/${REFDIR}/${REQSTEM}${_zenodo_ext}")

endfunction(octave_compare)


function(matlab_compare TESTNAME OUTDIR REFDIR REQSTEM)

add_test(NAME ${TESTNAME}
  COMMAND ${Matlab_MAIN_PROGRAM} -batch compare_all('${CMAKE_BINARY_DIR}/${OUTDIR}','${CMAKE_SOURCE_DIR}/${REFDIR}')
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 60
  REQUIRED_FILES "${CMAKE_BINARY_DIR}/${OUTDIR}/${REQSTEM}${OUTEXT};${CMAKE_SOURCE_DIR}/${REFDIR}/${REQSTEM}${_zenodo_ext}"
)
# Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endfunction(matlab_compare)
