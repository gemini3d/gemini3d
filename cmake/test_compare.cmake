
function(matlab_compare outdir refdir testname)

if(hdf5 OR netcdf)

add_test(NAME gemini:compare:${testname}:matlab
COMMAND ${Matlab_MAIN_PROGRAM} -batch "gemini3d.compare('${outdir}', '${refdir}')"
WORKING_DIRECTORY ${MATGEMINI_DIR}
)

set_tests_properties(gemini:compare:${testname}:matlab PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED "hdf5:${testname};netcdf:${testname}"
REQUIRED_FILES ${outdir}/inputs/config.nml)

endif()

endfunction(matlab_compare)


function(python_compare outdir refdir testname)

if(hdf5)
add_test(NAME gemini:compare:hdf5:${testname}:python
COMMAND ${Python3_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format h5
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

set_tests_properties(gemini:compare:hdf5:${testname}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED hdf5:${testname}
REQUIRED_FILES ${outdir}/inputs/config.nml)
endif(hdf5)

if(netcdf)
add_test(NAME gemini:compare:netcdf:${testname}:python
COMMAND ${Python3_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format nc)

set_tests_properties(gemini:compare:netcdf:${testname}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED netcdf:${testname}
REQUIRED_FILES ${outdir}/inputs/config.nml)
endif(netcdf)

endfunction(python_compare)


function(fortran_compare outdir refdir testname)

if(hdf5)

  add_test(NAME gemini:compare:hdf5:${testname}
  COMMAND $<TARGET_FILE:gemini3d.compare> ${outdir} ${refdir})

  set_tests_properties(gemini:compare:hdf5:${testname} PROPERTIES
  TIMEOUT 60
  FIXTURES_REQUIRED hdf5:${testname}
  REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml")

endif(hdf5)

if(netcdf)

  add_test(NAME gemini:compare:netcdf:${testname}
  COMMAND $<TARGET_FILE:gemini3d.compare> ${outdir} ${refdir})

  set_tests_properties(gemini:compare:netcdf:${testname} PROPERTIES
  TIMEOUT 60
  FIXTURES_REQUIRED netcdf:${testname}
  REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml")

endif(netcdf)

endfunction(fortran_compare)


function(compare_gemini_output testname)

set(outdir ${PROJECT_BINARY_DIR}/test${testname})
set(refdir ${PROJECT_SOURCE_DIR}/tests/data/test${testname})

if(MATGEMINI_DIR)
  matlab_compare(${outdir} ${refdir} ${testname})
endif(MATGEMINI_DIR)

if(PYGEMINI_DIR)
  python_compare(${outdir} ${refdir} ${testname})
endif(PYGEMINI_DIR)

# we always test this as it's sub-second runtime
fortran_compare(${outdir} ${refdir} ${testname})

endfunction(compare_gemini_output)
