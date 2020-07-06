include(${CMAKE_CURRENT_LIST_DIR}/matlab.cmake)


function(matlab_compare outdir refdir testname)

if(NOT matlab_ok)
  return()
endif()


add_test(NAME gemini:compare:hdf5:${testname}:matlab
  COMMAND ${Matlab_MAIN_PROGRAM} -batch "setup; compare_all('${outdir}', '${refdir}')"
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

set_tests_properties(gemini:compare:hdf5:${testname}:matlab PROPERTIES
TIMEOUT 120
DEPENDS gemini:hdf5:${testname}
REQUIRED_FILES ${outdir}/inputs/config.nml
SKIP_RETURN_CODE 77)

endfunction(matlab_compare)


function(python_compare outdir refdir testname)

if(python_ok AND hdf5)
  add_test(NAME gemini:compare:hdf5:${testname}:python
  COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gemini_compare.py ${outdir} ${refdir} -file_format h5
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

  set_tests_properties(gemini:compare:hdf5:${testname}:python PROPERTIES
  TIMEOUT 60
  DEPENDS gemini:hdf5:${testname}
  REQUIRED_FILES ${outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()

if(python_ok AND netcdf)
  add_test(NAME gemini:compare:netcdf:${testname}:python
  COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/gemini_compare.py ${outdir} ${refdir} -file_format nc)

  set_tests_properties(gemini:compare:netcdf:${testname}:python PROPERTIES
  TIMEOUT 60
  DEPENDS gemini:netcdf:${testname}
  REQUIRED_FILES ${outdir}/inputs/config.nml
  SKIP_RETURN_CODE 77)
endif()

endfunction(python_compare)


function(compare_gemini_output testname)
  set(outdir ${PROJECT_BINARY_DIR}/test${testname})
  set(refdir ${PROJECT_SOURCE_DIR}/tests/data/test${testname})

  matlab_compare(${outdir} ${refdir} ${testname})
  python_compare(${outdir} ${refdir} ${testname})
endfunction(compare_gemini_output)
