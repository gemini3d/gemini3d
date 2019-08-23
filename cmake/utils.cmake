include(${CMAKE_CURRENT_LIST_DIR}/maxfactor.cmake)

function(download_testfiles HASH REFNUM REFNAME ROOT)
# FetchContent is too aggressive, it deletes the output directory before extracting--could delete wrong directory causing data loss.
# So use this more specific but safe method.

set(ARCHIVE ${REFNAME}.zip)
set(URL https://zenodo.org/record/${REFNUM}/files/${REFNAME}.zip?download=1)

# --- ensure reference data is available for self-test
if(NOT EXISTS ${ROOT}/${REFNAME})

  if(EXISTS ${ROOT}/${ARCHIVE})
    file(MD5 ${ROOT}/${ARCHIVE} FHASH)
  endif()

  if(NOT EXISTS ${ROOT}/${ARCHIVE} OR NOT FHASH STREQUAL HASH)
    file(DOWNLOAD ${URL} ${ROOT}/${ARCHIVE}
      SHOW_PROGRESS)
	  # don't check hash during download in case of failure, so it doesn't stop using rest of program.

	file(MD5 ${ROOT}/${ARCHIVE} FHASH)

	# have user manually download
	if(NOT FHASH STREQUAL HASH)
	  message(WARNING "Self-tests may not work. Reference data did not match hash ${HASH} \n Download ${URL}  to  ${ARCHIVE} \n")
	endif()
  endif()

  # extract zip file (os-agnostic)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${ARCHIVE}
                  WORKING_DIRECTORY ${ROOT})
endif()

endfunction(download_testfiles)


function(num_mpi_processes REFDIR)

if(DEFINED NP)  # if the user manually sets NP at command line, don't override the user.
  return()
endif()

if(NOT WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
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

  set(NP ${MAXFACTOR} PARENT_SCOPE)
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 4)
  set(NP 4 PARENT_SCOPE)
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 2)
  set(NP 2 PARENT_SCOPE)
else()
  set(NP 1 PARENT_SCOPE)
endif()

endfunction(num_mpi_processes)


function(setup_gemini_test TESTNAME EXE TESTDIR REFDIR TIMEOUT)

num_mpi_processes(${REFDIR})

if(NP EQUAL 1)
  message(WARNING "Gemini with less than two MPI processes gives incorrect results on some MPI setups. We urge confirming results on a system with at least two MPI processes.")
endif()

set(TESTNAME ${TESTNAME}-NP${NP})  # for convenience, name with number of processes since this is important for debugging MPI

add_test(NAME ${TESTNAME}
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:${EXE}> ${CMAKE_CURRENT_SOURCE_DIR}/initialize/${TESTDIR}/config.ini ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  REQUIRED_FILES ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini
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
