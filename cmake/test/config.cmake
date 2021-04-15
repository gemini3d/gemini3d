include(${CMAKE_CURRENT_LIST_DIR}/compare.cmake)


function(setup_gemini_test name TIMEOUT)

# --- setup test
set(outdir ${PROJECT_BINARY_DIR}/test${name})
set(refroot ${PROJECT_SOURCE_DIR}/test_data)
set(refdir ${refroot}/test${name})

add_test(NAME gemini:${name}:setup
  COMMAND ${CMAKE_COMMAND} -Dname=${name} -Doutdir:PATH=${outdir} -Drefroot:PATH=${refroot} -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/download.cmake)
set_tests_properties(gemini:${name}:setup PROPERTIES
  FIXTURES_SETUP ${name}:setup
  FIXTURES_REQUIRED gemini_exe_fxt
  LABELS setup
  TIMEOUT 180)

# construct command
set(test_cmd $<TARGET_FILE:gemini3d.run> ${outdir} -gemexe $<TARGET_FILE:gemini.bin>)

if(mpi)
  list(APPEND test_cmd -mpiexec ${MPIEXEC_EXECUTABLE})

  if(NOT HWLOC_FOUND)
    # hwloc is probably the most accurate way to determine CPU count--fallback to CMake count.
    list(APPEND test_cmd -n ${Ncpu})
  endif()
endif()



if(hdf5)

add_test(NAME gemini:hdf5:${name}:dryrun
  COMMAND ${test_cmd} -dryrun)
# we prefer default WorkingDirectory of PROJECT_BINARY_DIR to make MSIS 2.0 msis20.parm use simpler
# otherwise, we have to generate source for msis_interface.f90

set_tests_properties(gemini:hdf5:${name}:dryrun PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED ${name}:setup
  FIXTURES_SETUP hdf5:${name}:dryrun
  REQUIRED_FILES ${outdir}/inputs/config.nml
  LABELS core)


add_test(NAME gemini:hdf5:${name} COMMAND ${test_cmd})

set_tests_properties(gemini:hdf5:${name} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED hdf5:${name}:dryrun
  FIXTURES_SETUP hdf5:${name}
  LABELS core)

endif(hdf5)


if(netcdf)
add_test(NAME gemini:netcdf:${name}:dryrun
  COMMAND ${test_cmd} -out_format nc -dryrun)

set_tests_properties(gemini:netcdf:${name}:dryrun PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED "mumps_fxt;${name}:setup"
  FIXTURES_SETUP netcdf:${name}:dryrun
  REQUIRED_FILES ${outdir}/inputs/config.nml
  LABELS core)

add_test(NAME gemini:netcdf:${name}
  COMMAND ${test_cmd} -out_format nc)

set_tests_properties(gemini:netcdf:${name} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED netcdf:${name}:dryrun
  FIXTURES_SETUP netcdf:${name}
  LABELS core)
endif(netcdf)

compare_gemini_output(${name} ${outdir} ${refdir})

endfunction(setup_gemini_test)
