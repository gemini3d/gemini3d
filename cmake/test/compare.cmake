function(python_compare outdir refdir name)

if(hdf5)

add_test(NAME gemini:compare:hdf5:${name}:python
COMMAND ${Python3_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format h5)

set_tests_properties(gemini:compare:hdf5:${name}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED hdf5:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS "compare;python")

endif(hdf5)

if(netcdf)

add_test(NAME gemini:compare:netcdf:${name}:python
COMMAND ${Python3_EXECUTABLE} -m gemini3d.compare ${outdir} ${refdir} -file_format nc)

set_tests_properties(gemini:compare:netcdf:${name}:python PROPERTIES
TIMEOUT 120
FIXTURES_REQUIRED netcdf:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS "compare;python")
endif(netcdf)

endfunction(python_compare)


function(fortran_compare outdir refdir name)

if(hdf5)

add_test(NAME gemini:compare:hdf5:${name}
COMMAND $<TARGET_FILE:gemini3d.compare> ${outdir} ${refdir})

set_tests_properties(gemini:compare:hdf5:${name} PROPERTIES
TIMEOUT 60
FIXTURES_REQUIRED hdf5:${name}:run_fxt
RESOURCE_LOCK cpu_mpi
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS compare)

# resource_lock compare for Windows, which can take 100x longer when run
# at same time with non-dependent sim runs.
# it's not a problem to run multiple compare at once, but it is a problem
# to run gemini3d.compare at same time as gemini.bin, even on different sims

endif(hdf5)

if(netcdf)

add_test(NAME gemini:compare:netcdf:${name}
COMMAND $<TARGET_FILE:gemini3d.compare> ${outdir} ${refdir})

set_tests_properties(gemini:compare:netcdf:${name} PROPERTIES
TIMEOUT 60
FIXTURES_REQUIRED netcdf:${name}:run_fxt
REQUIRED_FILES "${outdir}/inputs/config.nml;${refdir}/inputs/config.nml"
LABELS compare)

endif(netcdf)

endfunction(fortran_compare)


function(compare_gemini_output name outdir refdir)

if(PYGEMINI_DIR)
  python_compare(${outdir} ${refdir} ${name})
endif(PYGEMINI_DIR)

if(TARGET gemini3d.compare)
  fortran_compare(${outdir} ${refdir} ${name})
endif()

endfunction(compare_gemini_output)
