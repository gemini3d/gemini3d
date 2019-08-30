function(num_mpi_processes REFDIR)

if(DEFINED NP)  # if the user manually sets NP at command line, don't override the user.
  return()
endif()

if(NOT WIN32)
  set(SIZEFN ${REFDIR}/inputs/simsize.dat)

  # file(READ ${SIZEFN} hex OFFSET 8 LIMIT 4 HEX) # no, this leaves trailing zeros
  execute_process(
    COMMAND od -j4 -N9 -t d4 ${SIZEFN}
    OUTPUT_VARIABLE raw
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  separate_arguments(x2x3 UNIX_COMMAND ${raw})
  # message(STATUS "x2x3 ${x2x3}")
  list(GET x2x3 1 x2)
  list(GET x2x3 2 x3)

  if(x3 EQUAL 1)  # 2D sim
    set(x3 ${x2})
  endif()

  #message(STATUS "x2 ${x2} x3 ${x3}")
  math(EXPR halfX3 "${x3} / 2")
  #message("halfx3 ${halfX3}")

  maxfactor(${halfX3} ${MPIEXEC_MAX_NUMPROCS})

  set(NP ${MAXFACTOR})
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 4)
  set(NP 4)
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 2)
  set(NP 2)
else()  # Gemini needs at least 2 MPI images, even if only single CPU core available.
  set(NP 2)
endif()

# Gemini needs at least 2 MPI images, even if only single CPU core available.
if(NP LESS 2)
  set(NP 2)
endif()

set(NP ${NP} PARENT_SCOPE)

endfunction(num_mpi_processes)