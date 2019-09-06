include(${CMAKE_CURRENT_LIST_DIR}/maxfactor.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/num_mpi.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/download.cmake)

function(setup_gemini_test TESTNAME EXE TESTDIR REFDIR TIMEOUT)

num_mpi_processes(${REFDIR})

if(NP LESS 2)
  message(FATAL_ERROR "Gemini with less than two MPI processes will fail. Must use at least 2 MPI images, even on single core CPU.")
endif()

set(TESTNAME ${TESTNAME}-NP${NP})  # for convenience, name with number of processes since this is important for debugging MPI

add_test(NAME ${TESTNAME}
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${EXE}> ${CMAKE_CURRENT_SOURCE_DIR}/initialize/${TESTDIR}/config.nml ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  REQUIRED_FILES ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.nml
  FIXTURES_REQUIRED "MPIMUMPS;IOfmt"
  RUN_SERIAL true
)

endfunction(setup_gemini_test)


function(compare_gemini_output TESTNAME OUTDIR REFDIR REQFILE)
# This sets up the Compare* tests

include(${PROJECT_SOURCE_DIR}/cmake/compare.cmake)

#--- Python
if(PythonOK)
  python_compare(${TESTNAME}_Python ${OUTDIR} ${REFDIR} ${REQFILE})
endif(PythonOK)

#--- Octave
if(OctaveOK)
  octave_compare(${TESTNAME}_Octave ${OUTDIR} ${REFDIR} ${REQFILE})
endif()

#--- Matlab
if(MatlabOK)
  matlab_compare(${TESTNAME}_Matlab ${TESTDIR} ${REFDIR} ${REQFILE})
endif()

endfunction(compare_gemini_output)
