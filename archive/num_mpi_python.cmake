function(num_mpi_processes REFDIR)

set(_pyok)

if(PythonOK)
  # do not quote COMMAND line in general
  set(_cmd ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/script_utils/get_cpu_count.py ${REFDIR}/inputs)

  execute_process(COMMAND ${_cmd}
    RESULT_VARIABLE _pyok
    OUTPUT_VARIABLE NP
    TIMEOUT 15)
endif()

if(NOT _pyok EQUAL 0 OR NP LESS 1)
  message(STATUS "${_cmd} could not detect MPI count, falling back to 1 (tests will run slowly)")
  set(NP 1)
endif()

set(NP ${NP} PARENT_SCOPE)

endfunction(num_mpi_processes)

# test standalone:
#   cmake -P num_mpi_processes(tests/data/zenodo2d_fang)