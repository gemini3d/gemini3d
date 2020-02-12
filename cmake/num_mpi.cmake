function(num_mpi_processes REFDIR)

if(PythonOK)
# do not quote COMMAND line in general
execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/script_utils/get_cpu_count.py ${REFDIR}/inputs
  OUTPUT_VARIABLE NP
  TIMEOUT 15)
else()
  # CMake's CPU count is not so reliable
  include(ProcessorCount)
  ProcessorCount(NP)
endif()

if(NP LESS 1)
  message(STATUS "could not detect number of CPUs, falling back to 1")
  set(NP 1)
endif()

set(NP ${NP} PARENT_SCOPE)

endfunction(num_mpi_processes)

# for testing standalone cmake -P
# num_mpi_processes(tests/data/zenodo2d_fang)