include(${CMAKE_CURRENT_LIST_DIR}/num_mpi.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/download.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/compare.cmake)

function(setup_gemini_test TESTNAME EXE TESTDIR REFDIR TIMEOUT)

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR})
  message(WARNING "directory ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR} not found.")
  return()
endif()

# sets NP
if(NOT NP)
  num_mpi_processes(${REFDIR})
endif()

set(_name ${TESTNAME})

add_test(NAME ${_name}
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${EXE}> ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.nml ${CMAKE_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${_name} PROPERTIES
  TIMEOUT ${TIMEOUT}
  SKIP_RETURN_CODE 77
  FIXTURES_REQUIRED "MPIMUMPS;IOfmt"
  RUN_SERIAL true
)
if(WIN32 AND HDF5_ROOT)  # for Windows ifort dll
  set_tests_properties(${_name} PROPERTIES ENVIRONMENT "PATH=${HDF5_ROOT}/bin;$ENV{PATH}")
endif()

endfunction(setup_gemini_test)
