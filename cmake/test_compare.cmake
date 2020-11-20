
function(matlab_compare outdir refdir testname)

if(hdf5 OR netcdf)

add_test(NAME gemini:compare:${testname}:matlab
COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${CMAKE_CURRENT_SOURCE_DIR}/setup.m'); gemini3d.compare_all('${outdir}', '${refdir}')"
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
TIMEOUT 60
FIXTURES_REQUIRED hdf5:${testname}
REQUIRED_FILES ${outdir}/inputs/config.nml)
endif(hdf5)

if(netcdf)
add_test(NAME gemini:compare:netcdf:${testname}:python
COMMAND ${Python3_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format nc)

set_tests_properties(gemini:compare:netcdf:${testname}:python PROPERTIES
TIMEOUT 60
FIXTURES_REQUIRED netcdf:${testname}
REQUIRED_FILES ${outdir}/inputs/config.nml)
endif(netcdf)

endfunction(python_compare)


function(compare_gemini_output testname)

set(outdir ${PROJECT_BINARY_DIR}/test${testname})
set(refdir ${PROJECT_SOURCE_DIR}/tests/data/test${testname})

if(matlab)
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if(Matlab_FOUND)
matlab_compare(${outdir} ${refdir} ${testname})
endif(Matlab_FOUND)
endif(matlab)

if(python_ok)
python_compare(${outdir} ${refdir} ${testname})
endif()

endfunction(compare_gemini_output)
