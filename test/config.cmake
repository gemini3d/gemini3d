# Audit modification 2026-09-16: correct MPI test worker count and fail if required reference tests cannot be registered
set(gemini3d_test_run_root "${PROJECT_BINARY_DIR}/test_runs" CACHE PATH "Isolated native test output directory")
if(NOT IS_ABSOLUTE "${gemini3d_test_run_root}")
  message(FATAL_ERROR "gemini3d_test_run_root must be absolute")
endif()
function(gcd a b out_var)
  if(a LESS 1 OR b LESS 1)
    message(FATAL_ERROR "gcd: positive integers only")
  endif()

  set(x ${a})
  set(y ${b})
  math(EXPR z "${x} % ${y}")

  while(NOT z EQUAL 0)
    set(x ${y})
    set(y ${z})
    math(EXPR z "${x} % ${y}")
  endwhile()

  set(${out_var} ${y} PARENT_SCOPE)
endfunction()


function(max_gcd L M out_var)
## largest divisor of L that is <= M
if(M LESS 1)
  message(FATAL_ERROR "max_gcd: CPU count must be at least one")
endif()
set(i ${M})
while(i GREATER_EQUAL 1)
  math(EXPR r "${L} % ${i}")
  if(r EQUAL 0)
    set(${out_var} ${i} PARENT_SCOPE)
    return()
  endif()
  math(EXPR i "${i} - 1")
endwhile()
endfunction()


function(max_gcd2 lx2 lx3 M out_var)
## find (d2,d3) with d2|lx2, d3|lx3, d2*d3<=M that maximises d2*d3,
## return the worker count; process-grid factor selection happens in the solver
if(M LESS 1)
  message(FATAL_ERROR "max_gcd2: CPU count must be at least one")
endif()

if(lx2 LESS 1 OR lx3 LESS 1)
  message(FATAL_ERROR "max_gcd2: grid dimensions must be positive")
endif()
set(best 1)
math(EXPR limit "${M}")
if(lx2 LESS limit)
  set(limit ${lx2})
endif()
foreach(d2 RANGE 1 ${limit})
  math(EXPR remainder "${lx2} % ${d2}")
  if(NOT remainder EQUAL 0)
    continue()
  endif()
  math(EXPR budget "${M} / ${d2}")
  max_gcd(${lx3} ${budget} d3)
  math(EXPR candidate "${d2} * ${d3}")
  if(candidate GREATER best)
    set(best ${candidate})
  endif()
endforeach()
set(${out_var} ${best} PARENT_SCOPE)

endfunction()


function(max_mpi lx2 lx3 max_cpu out_var)
## return max useful MPI worker count for an lx2 x lx3 simulation grid
math(EXPR _lx2h "${lx2} / 2")
math(EXPR _lx3h "${lx3} / 2")

if(lx3 EQUAL 1)
  max_gcd(${_lx2h} ${max_cpu} _n)
elseif(lx2 EQUAL 1)
  max_gcd(${_lx3h} ${max_cpu} _n)
else()
  max_gcd2(${_lx2h} ${_lx3h} ${max_cpu} _n)
endif()

set(${out_var} ${_n} PARENT_SCOPE)
endfunction()


function(setup_gemini_test name)

include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/url_name.cmake)

if(name MATCHES "_cpp$" AND NOT TARGET gemini_c.bin)
  return()
endif()

# --- setup test
set(out_dir ${gemini3d_test_run_root}/${name})
set(ref_root ${PROJECT_BINARY_DIR}/test_data/compare)
set(ref_dir ${ref_root}/${name})
set(arc_json_file ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/ref_data_snapshot.json)
if(NOT EXISTS "${arc_json_file}")
  message(FATAL_ERROR "Required pinned reference manifest is missing: ${arc_json_file}")
endif()

# --- compute proper number of MPI workers for this test

file(READ ${arc_json_file} _refj)

get_url_name(${name} url_name)

string(JSON Nlx ERROR_VARIABLE _err LENGTH "${_refj}" tests ${url_name} lx)
if(_err)
  message(FATAL_ERROR "test ${name} missing lx in ${arc_json_file}: ${_err}")
  return()
endif()
if(NOT Nlx EQUAL 3)
  message(FATAL_ERROR "test ${name} has lx=${Nlx} in ${arc_json_file}, expected 3")
  return()
endif()

string(JSON lx1 GET "${_refj}" tests ${url_name} lx 0)
string(JSON lx2 GET "${_refj}" tests ${url_name} lx 1)
string(JSON lx3 GET "${_refj}" tests ${url_name} lx 2)

max_mpi(${lx2} ${lx3} ${MPIEXEC_MAX_NUMPROCS} Nworker)

message(STATUS "test ${name}: lx=(${lx1},${lx2},${lx3}) Nworker=${Nworker}")

# --- define tests

add_test(NAME ${name}:download
COMMAND ${CMAKE_COMMAND}
  -Dname=${name}
  -Doutdir:PATH=${out_dir}
  -Drefroot:PATH=${ref_root}
  -Darc_json_file:FILEPATH=${arc_json_file}
  -Dreset_test_outputs=ON
  -Dtest_binary_dir:PATH=${gemini3d_test_run_root}
  -P ${CMAKE_CURRENT_LIST_DIR}/download.cmake
)
set_tests_properties(${name}:download PROPERTIES
FIXTURES_SETUP ${name}:download_fxt
RESOURCE_LOCK download_lock  # avoid anti-leeching transient failures
LABELS download
)

# --- gemini3d.run ---
set(test_cmd gemini3d.run ${out_dir} -mpiexec ${MPIEXEC_EXECUTABLE})
if(name MATCHES "_cpp$")
  list(APPEND test_cmd -exe $<TARGET_FILE:gemini_c.bin>)
else()
  list(APPEND test_cmd -exe $<TARGET_FILE:gemini.bin>)
endif()

add_test(NAME gemini_run:${name}:dryrun COMMAND ${test_cmd} -dryrun)
set_tests_properties(gemini_run:${name}:dryrun PROPERTIES
FIXTURES_REQUIRED "gemini_exe_fxt;${name}:download_fxt"
FIXTURES_SETUP ${name}:frontend_dryrun
RESOURCE_LOCK cpu_mpi
WORKING_DIRECTORY $<TARGET_FILE_DIR:gemini3d.run>
PROCESSORS ${Nworker}
)
hdf5_dll(gemini_run:${name}:dryrun)

# --- gemini.bin dryrun ---

if(name MATCHES "_cpp$")
  set(test_cmd gemini_c.bin)
else()
  set(test_cmd gemini.bin)
endif()

# $<TARGET_FILE:${test_cmd}> is essential when using direct command line in add_test

test_mpi_command(${Nworker} "$<TARGET_FILE_DIR:gemini.bin>" mpi_cmd)
add_test(NAME gemini:${name}:dryrun COMMAND ${mpi_cmd} $<TARGET_FILE:${test_cmd}> ${out_dir} -dryrun)
test_mpi_props(gemini:${name}:dryrun ${Nworker})
set_tests_properties(gemini:${name}:dryrun PROPERTIES
FIXTURES_SETUP ${name}:dryrun
FIXTURES_REQUIRED "gemini_exe_fxt;${name}:download_fxt;${name}:frontend_dryrun"
)
hdf5_dll(gemini:${name}:dryrun)

# --- gemini.bin run ---

test_mpi_command(${Nworker} "$<TARGET_FILE_DIR:gemini.bin>" mpi_cmd)
add_test(NAME gemini:${name} COMMAND ${mpi_cmd} $<TARGET_FILE:${test_cmd}> ${out_dir})
test_mpi_props(gemini:${name} ${Nworker})
set_tests_properties(gemini:${name} PROPERTIES
FIXTURES_REQUIRED ${name}:dryrun
FIXTURES_SETUP ${name}:run_fxt
)
hdf5_dll(gemini:${name})
# WORKING_DIRECTORY is needed for tests like HWM14 that need data files in binary directory.
set_tests_properties(gemini:${name}:dryrun gemini:${name} PROPERTIES
RESOURCE_LOCK cpu_mpi
REQUIRED_FILES ${out_dir}/inputs/config.nml
LABELS core
)
if(name MATCHES "_cpp$")
  set_property(TEST gemini:${name}:dryrun gemini:${name} PROPERTY LABELS "gemini3d:core")
endif()

compare_gemini_output(${name} ${out_dir} ${ref_dir})

endfunction()


function(setup_magcalc_test name)

set(out_dir ${gemini3d_test_run_root}/${name})

add_test(NAME magcalc:${name}:setup
COMMAND ${Python_EXECUTABLE} -m gemini3d.magcalc ${out_dir}
)
set_tests_properties(magcalc:${name}:setup PROPERTIES
FIXTURES_REQUIRED ${name}:run_fxt
FIXTURES_SETUP magcalc:${name}:setup
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)
hdf5_dll(magcalc:${name}:setup)

add_test(NAME magcalc:${name} COMMAND magcalc.run ${out_dir})
set_tests_properties(magcalc:${name} PROPERTIES
RESOURCE_LOCK cpu_mpi
FIXTURES_REQUIRED magcalc:${name}:setup
LABELS core
DISABLED $<NOT:$<BOOL:${PYGEMINI_DIR}>>
)
hdf5_dll(magcalc:${name})

endfunction(setup_magcalc_test)
