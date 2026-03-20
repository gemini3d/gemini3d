function(setup_gemini_test name)

if(name MATCHES "_cpp$" AND NOT TARGET gemini_c.bin)
  return()
endif()

# --- setup test
set(out_dir ${PROJECT_BINARY_DIR}/${name})
set(ref_root ${PROJECT_BINARY_DIR}/test_data/compare)
set(ref_dir ${ref_root}/${name})
set(arc_json_file ${PROJECT_BINARY_DIR}/ref_data.json)

add_test(NAME ${name}:download
COMMAND ${CMAKE_COMMAND}
  -Dname=${name}
  -Doutdir:PATH=${out_dir}
  -Drefroot:PATH=${ref_root}
  -Darc_json_file:FILEPATH=${arc_json_file}
  -P ${CMAKE_CURRENT_LIST_DIR}/download.cmake
)
set_tests_properties(${name}:download PROPERTIES
FIXTURES_SETUP ${name}:download_fxt
RESOURCE_LOCK download_lock  # avoid anti-leeching transient failures
LABELS download
)

# construct command
set(test_cmd gemini3d.run ${out_dir})
if(name MATCHES "_cpp$")
  list(APPEND test_cmd -exe $<TARGET_FILE:gemini_c.bin>)
else()
  list(APPEND test_cmd -exe $<TARGET_FILE:gemini.bin>)
endif()
list(APPEND test_cmd -mpiexec ${MPIEXEC_EXECUTABLE})

add_test(NAME gemini:${name}:dryrun
COMMAND ${test_cmd} -dryrun
)

set_tests_properties(gemini:${name}:dryrun PROPERTIES
FIXTURES_SETUP ${name}:dryrun
FIXTURES_REQUIRED "gemini_exe_fxt;${name}:download_fxt"
)


add_test(NAME gemini:${name} COMMAND ${test_cmd})

set_tests_properties(gemini:${name} PROPERTIES
FIXTURES_REQUIRED ${name}:dryrun
FIXTURES_SETUP ${name}:run_fxt
)

# WORKING_DIRECTORY is needed for tests like HWM14 that need data files in binary directory.
set_tests_properties(gemini:${name}:dryrun gemini:${name} PROPERTIES
RESOURCE_LOCK cpu_mpi
REQUIRED_FILES ${out_dir}/inputs/config.nml
LABELS core
WORKING_DIRECTORY $<TARGET_FILE_DIR:gemini.bin>
)
if(name MATCHES "_cpp$")
  set_property(TEST gemini:${name}:dryrun gemini:${name} PROPERTY LABELS "core;Cpp")
endif()
if(DEFINED mpi_tmpdir)
  set_property(TEST gemini:${name}:dryrun gemini:${name} PROPERTY ENVIRONMENT "TMPDIR=${mpi_tmpdir}")
endif()


compare_gemini_output(${name} ${out_dir} ${ref_dir})

endfunction(setup_gemini_test)


function(setup_magcalc_test name)

set(out_dir ${PROJECT_BINARY_DIR}/${name})

add_test(NAME magcalc:${name}:setup
COMMAND ${Python_EXECUTABLE} -m gemini3d.magcalc ${out_dir}
)
set_tests_properties(magcalc:${name}:setup PROPERTIES
FIXTURES_REQUIRED ${name}:run_fxt
FIXTURES_SETUP magcalc:${name}:setup
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)

add_test(NAME magcalc:${name} COMMAND magcalc.run ${out_dir})
set_tests_properties(magcalc:${name} PROPERTIES
RESOURCE_LOCK cpu_mpi
FIXTURES_REQUIRED magcalc:${name}:setup
LABELS core
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)

endfunction(setup_magcalc_test)
