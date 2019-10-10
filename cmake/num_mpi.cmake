function(num_mpi_processes REFDIR)

if(DEFINED NP)  # if the user manually sets NP at command line, don't override the user.
  return()
endif()

if(PythonOK)
  execute_process(COMMAND ${Python3_EXECUTABLE} script_utils/meson_cpu_count.py ${REFDIR}/inputs/simsize.dat
    ERROR_VARIABLE oops
    OUTPUT_VARIABLE NP
    RESULT_VARIABLE ret
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    TIMEOUT 15)
  if(NOT ret EQUAL 0)
    message(WARNING "${oops}")
    set(NP 2)
  endif()
else()
  set(NP 2)
endif(PythonOK)

# Gemini needs at least 2 MPI images, even if only single CPU core available.
if(NP LESS 2)
  set(NP 2)
endif()

set(NP ${NP} PARENT_SCOPE)

endfunction(num_mpi_processes)