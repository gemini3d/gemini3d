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

	# try 2nd time using curl
	if(NOT FHASH STREQUAL HASH)
	  execute_process(COMMAND curl -L --url ${URL} -o ${ARCHIVE} WORKING_DIRECTORY ${ROOT})
	endif()

	file(MD5 ${ROOT}/${ARCHIVE} FHASH)
	# have user manually download
	if(NOT FHASH STREQUAL HASH)
	  message(WARNING "Self-tests may not work. Reference data did not match hash ${HASH} \n Download ${URL}  to  ${ARCHIVE} \n")
	endif()
  endif()

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

  message(STATUS "Gemini test auto-setup with ${MAXFACTOR} MPI processes")

  set(NP ${MAXFACTOR} PARENT_SCOPE)
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 4)
  set(NP 4 PARENT_SCOPE)
elseif(${MPIEXEC_MAX_NUMPROCS} GREATER_EQUAL 2)
  set(NP 2 PARENT_SCOPE)
else()
  set(NP 1 PARENT_SCOPE)
endif()

endfunction(num_mpi_processes)


function(check_octave_source_runs code)

if(NOT Octave_EXECUTABLE)
  set(ok false)
else()
  execute_process(COMMAND ${Octave_EXECUTABLE} --eval ${code}
    ERROR_QUIET OUTPUT_QUIET
    RESULT_VARIABLE ret
    TIMEOUT 5)
  if(ret EQUAL 0)
    set(ok true)
  else()
    set(ok false)
  endif()
endif()

set(OctaveOK ${ok} CACHE BOOL "GNU Octave is sufficiently new to run self-tests")

endfunction(check_octave_source_runs)


function(check_matlab_source_runs code)

if(NOT Matlab_MAIN_PROGRAM)
  set(ok false)
else()
  execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch ${code}
    ERROR_QUIET OUTPUT_QUIET
    RESULT_VARIABLE ret
    TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes
  if(ret EQUAL 0)
    set(ok true)
  else()
    set(ok false)
  endif()
endif()

set(MatlabOK ${ok} CACHE BOOL "Matlab is sufficiently new to run self-tests")

endfunction(check_matlab_source_runs)


function(setup_gemini_test TESTNAME TESTDIR REFDIR TIMEOUT)

num_mpi_processes(${REFDIR})

if(NP EQUAL 1)
  message(WARNING "Gemini with less than two MPI processes gives incorrect results on some MPI setups. We urge confirming results on a system with at least two MPI processes.")
endif()

set(TESTNAME ${TESTNAME}-NP${NP})  # for convenience, name with number of processes since this is important for debugging MPI

message(STATUS "Test ${TESTNAME} uses ${NP} MPI processes")

add_test(NAME ${TESTNAME}
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} $<TARGET_FILE:gemini.bin> ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT ${TIMEOUT}
  REQUIRED_FILES ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini
  FIXTURES_REQUIRED "MPIMUMPS;IOfmt"
  RUN_SERIAL true
)

endfunction(setup_gemini_test)


function(octave_compare TESTNAME OUTDIR REFDIR REQFILE)

add_test(NAME ${TESTNAME}
  COMMAND ${Octave_EXECUTABLE} --eval "exit(compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}'))"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 30
  REQUIRED_FILES "${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/${REQFILE};${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}/${REQFILE}"
)

endfunction(octave_compare)


function(matlab_compare TESTNAME OUTDIR REFDIR REQFILE)

add_test(NAME ${TESTNAME}
  COMMAND ${Matlab_MAIN_PROGRAM} -batch compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}')
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES
  TIMEOUT 60
  REQUIRED_FILES "${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/${REQFILE};${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}/${REQFILE}"
)
# Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endfunction(matlab_compare)


function(compare_gemini_output TESTNAME TESTDIR REFDIR REQFILE)

if(NOT DEFINED OctaveOK)
  if(WIN32)
    if(NOT Octave_ROOT)
      set(Octave_ROOT $ENV{HOMEDRIVE}/Octave)
    endif()
    file(GLOB _octpath "${Octave_ROOT}/Octave*/mingw64/bin/")
  endif()

  find_program(Octave_EXECUTABLE
    NAMES octave-cli octave
    DOC "GNU Octave >= 4.0"
    PATHS ${Octave_ROOT}
    HINTS ${_octpath}
    PATH_SUFFIXES bin)

  # https://octave.sourceforge.io/octave/function/exist.html
  check_octave_source_runs("exit(assert(exist('validateattributes', 'file')==2))")
endif()

if(OctaveOK)
  message(STATUS "Using GNU Octave: ${Octave_EXECUTABLE}")
  octave_compare(${TESTNAME} ${TESTDIR} ${REFDIR} ${REQFILE})
  return()
endif()

# --- try Matlab if Octave wasn't available
find_package(Matlab COMPONENTS MAIN_PROGRAM)

if (Matlab_FOUND)
  message(STATUS "Using Matlab: ${Matlab_MAIN_PROGRAM}")
  matlab_compare(Matlab${TESTNAME} ${TESTDIR} ${REFDIR} ${REQFILE})
else()
  message(WARNING "Neither Matlab or Octave was found, cannot run full self-test")
endif()

endfunction(compare_gemini_output)
