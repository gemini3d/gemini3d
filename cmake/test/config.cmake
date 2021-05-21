include(${CMAKE_CURRENT_LIST_DIR}/compare.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/system_meta.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/RefPackage.cmake)


function(setup_gemini_test name TIMEOUT)

# --- setup test
cmake_path(APPEND out_dir ${PROJECT_BINARY_DIR} ${name})
cmake_path(APPEND ref_root ${PROJECT_SOURCE_DIR} test_data)
cmake_path(APPEND ref_dir ${ref_root} ${name})

add_test(NAME ${name}:download
  COMMAND ${CMAKE_COMMAND} -Dname=${name} -Doutdir:PATH=${out_dir} -Drefroot:PATH=${ref_root} -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/download.cmake)
set_tests_properties(${name}:download PROPERTIES
  FIXTURES_SETUP ${name}:download_fxt
  LABELS download
  TIMEOUT 180)

# construct command
set(test_cmd $<TARGET_FILE:gemini3d.run> ${out_dir} -exe $<TARGET_FILE:gemini.bin>)

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
  FIXTURES_REQUIRED "gemini_exe_fxt;${name}:download_fxt"
  FIXTURES_SETUP hdf5:${name}:dryrun
  REQUIRED_FILES ${out_dir}/inputs/config.nml
  LABELS core)


add_test(NAME gemini:hdf5:${name} COMMAND ${test_cmd})

set_tests_properties(gemini:hdf5:${name} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED hdf5:${name}:dryrun
  FIXTURES_SETUP hdf5:${name}:run_fxt
  LABELS core)

endif(hdf5)


if(netcdf)
add_test(NAME gemini:netcdf:${name}:dryrun
  COMMAND ${test_cmd} -out_format nc -dryrun)

set_tests_properties(gemini:netcdf:${name}:dryrun PROPERTIES
  TIMEOUT 60
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED "gemini_exe_fxt;${name}:download_fxt"
  FIXTURES_SETUP netcdf:${name}:dryrun
  REQUIRED_FILES ${out_dir}/inputs/config.nml
  LABELS core)

add_test(NAME gemini:netcdf:${name}
  COMMAND ${test_cmd} -out_format nc)

set_tests_properties(gemini:netcdf:${name} PROPERTIES
  TIMEOUT ${TIMEOUT}
  RESOURCE_LOCK cpu_mpi
  FIXTURES_REQUIRED netcdf:${name}:dryrun
  FIXTURES_SETUP netcdf:${name}:run_fxt
  LABELS core)
endif(netcdf)

if(package)
  ref_package(${out_dir} ${ref_json_file} ${name})
else()
  compare_gemini_output(${name} ${out_dir} ${ref_dir})
endif()

endfunction(setup_gemini_test)
