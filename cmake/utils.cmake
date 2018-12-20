
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
      SHOW_PROGRESS
		  EXPECTED_HASH MD5=${HASH})
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${ARCHIVE}
                  WORKING_DIRECTORY ${ROOT})
endif()

endfunction()


function(check_ram)
# Cygwin does not correctly report memory or CPU resources to CMake_Host_system_info

if(CYGWIN)
  include(ProcessorCount)
  ProcessorCount(NCORES)
  math(EXPR NCORES "${NCORES}/2")
  set(NTEST ${NCORES} PARENT_SCOPE)
  return()
endif()

cmake_host_system_information(RESULT PHYSRAM QUERY AVAILABLE_PHYSICAL_MEMORY)

if(${PHYSRAM} LESS 1000)
  set(NTEST 1 PARENT_SCOPE)
else()
  set(NTEST ${MPIEXEC_MAX_NUMPROCS} PARENT_SCOPE)
endif()

endfunction(check_ram)


function(check_octave_source_runs code)

# imported target doesn't work with CMake 3.13
execute_process(COMMAND ${Octave_EXECUTABLE} --eval ${code}
  ERROR_QUIET
  RESULT_VARIABLE ok
  TIMEOUT 5)

set(OctaveOK ${ok} PARENT_SCOPE)

endfunction(check_octave_source_runs)


function(check_matlab_source_runs code)

execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -nojvm -r ${code}
  ERROR_QUIET OUTPUT_QUIET
  RESULT_VARIABLE ok
  TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes

set(MatlabOK ${ok} PARENT_SCOPE)

endfunction(check_matlab_source_runs)


function(run_gemini_test TESTNAME TESTDIR TIMEOUT)

check_ram()

add_test(NAME ${TESTNAME}
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NTEST} ${CMAKE_SOURCE_DIR}/gemini ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${TESTNAME} PROPERTIES 
  TIMEOUT ${TIMEOUT}
  REQUIRED_FILES ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini
)

endfunction()


function(octave_compare TESTNAME OUTDIR REFDIR REQFILE)

add_test(NAME ${TESTNAME}
  COMMAND Octave::Interpreter --eval "exit(compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}'))"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES 
  TIMEOUT 30
  REQUIRED_FILES "${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/${REQFILE};${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}/${REQFILE}"
)

endfunction()


function(matlab_compare TESTNAME OUTDIR REFDIR REQFILE)

add_test(NAME ${TESTNAME}
  COMMAND ${Matlab_MAIN_PROGRAM} -nojvm -r "exit(compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}'))"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES 
  TIMEOUT 60
  REQUIRED_FILES "${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}/${REQFILE};${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}/${REQFILE}"
) 
# Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endfunction()


function(compare_gemini_output TESTNAME TESTDIR REFDIR REQFILE)

if(NOT DEFINED OctaveOK)
  find_package(Octave COMPONENTS Interpreter)
  check_octave_source_runs("exit(exist('validateattributes'))")
endif()

if(OctaveOK)
  octave_compare(${TESTNAME} ${TESTDIR} ${REFDIR} ${REQFILE})
else()
  if(NOT DEFINED MatlabOK)
    find_package(Matlab QUIET COMPONENTS MAIN_PROGRAM)
    check_matlab_source_runs("exit(exist('validateattributes'))")
  endif()

  if (MatlabOK)
    matlab_compare(Matlab${TESTNAME} ${TESTDIR} ${REFDIR} ${REQFILE})
  else()
    message(WARNING "Neither Matlab or Octave was found, cannot run full self-test")
  endif()
endif()

endfunction()
