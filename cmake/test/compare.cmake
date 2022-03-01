
function(matlab_compare outdir refdir name)

if(hdf5 OR netcdf)

add_test(NAME gemini:compare:${name}:matlab
COMMAND ${Matlab_MAIN_PROGRAM} -batch "gemini3d.compare('${outdir}', '${refdir}')"
WORKING_DIRECTORY ${matgemini_SOURCE_DIR}
)

set_tests_properties(gemini:compare:${name}:matlab PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED "hdf5:${name}:run_fxt;netcdf:${name}:run_fxt"
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
ENVIRONMENT "MATLABPATH=${MATLABPATH}"
LABELS "compare;matlab"
)

endif()

endfunction(matlab_compare)


function(python_compare outdir refdir name)

add_test(NAME gemini:compare:hdf5:${name}:python
COMMAND ${Python_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format h5)

set_tests_properties(gemini:compare:hdf5:${name}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED hdf5:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS "compare;python"
DISABLED $<OR:$<NOT:$<BOOL:${PYGEMINI_DIR}>>,$<NOT:$<BOOL:${hdf5}>>>
)

if(netcdf)

add_test(NAME gemini:compare:netcdf:${name}:python
COMMAND ${Python_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format nc)

set_tests_properties(gemini:compare:netcdf:${name}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED netcdf:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS "compare;python"
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)

endif(netcdf)

endfunction(python_compare)


function(fortran_compare outdir refdir name)

add_test(NAME gemini:compare:hdf5:${name}
COMMAND gemini3d.compare ${outdir} ${refdir})

set_tests_properties(gemini:compare:hdf5:${name} PROPERTIES
TIMEOUT 60
FIXTURES_REQUIRED hdf5:${name}:run_fxt
RESOURCE_LOCK $<$<BOOL:${WIN32}>:cpu_mpi>
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS compare
DISABLED $<OR:$<NOT:$<TARGET_EXISTS:gemini3d.compare>>,$<NOT:$<BOOL:${hdf5}>>>
)
dll_test_path(h5fortran::h5fortran gemini:compare:hdf5:${name})

# resource_lock compare for Windows, which can take 100x longer when run
# at same time with non-dependent sim runs.
# it's not a problem to run multiple compare at once, but it is a problem
# to run gemini3d.compare at same time as gemini.bin, even on different sims

if(netcdf)

add_test(NAME gemini:compare:netcdf:${name}
COMMAND gemini3d.compare ${outdir} ${refdir})

set_tests_properties(gemini:compare:netcdf:${name} PROPERTIES
TIMEOUT 60
FIXTURES_REQUIRED netcdf:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS compare
)

endif(netcdf)

endfunction(fortran_compare)


function(compare_gemini_output name outdir refdir)

if(matlab)
  matlab_compare(${outdir} ${refdir} ${name})
endif()

if(python)
  python_compare(${outdir} ${refdir} ${name})
endif()

fortran_compare(${outdir} ${refdir} ${name})

endfunction(compare_gemini_output)
