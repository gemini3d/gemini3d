
function(download_testfiles HASH REFNUM REFNAME ROOT)

set(ARCHIVE ${REFNAME}.zip)
set(URL https://zenodo.org/record/${REFNUM}/files/${REFNAME}.zip?download=1)

# --- ensure 2D reference data is available for self-test
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


function(check_ram NTEST)

cmake_host_system_information(RESULT PHYSRAM QUERY AVAILABLE_PHYSICAL_MEMORY)
if(${PHYSRAM} LESS 1000)
  set(NTEST 1 PARENT_SCOPE)
else()
  set(NTEST 4 PARENT_SCOPE)
endif()

endfunction()


function(run_gemini_test TESTNAME TESTDIR TIMEOUT)

check_ram(NTEST)

add_test(NAME ${TESTNAME}
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NTEST} ${CMAKE_SOURCE_DIR}/gemini ${CMAKE_SOURCE_DIR}/initialize/${TESTDIR}/config.ini ${CMAKE_CURRENT_BINARY_DIR}/${TESTDIR}
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${TIMEOUT})

endfunction()


function(octave_compare TESTNAME OUTDIR REFDIR)

find_package(Octave COMPONENTS Interpreter)

if(Octave_Interpreter_FOUND)

  add_test(NAME ${TESTNAME}
           COMMAND Octave::Interpreter -q --eval "compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}')"
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 30)

endif()

endfunction()


function(matlab_compare TESTNAME OUTDIR REFDIR)

find_package(Matlab QUIET COMPONENTS MAIN_PROGRAM)

if (Matlab_MAIN_PROGRAM_FOUND)
add_test(NAME ${TESTNAME}
         COMMAND ${Matlab_MAIN_PROGRAM} -nojvm -r "exit(compare_all('${CMAKE_CURRENT_BINARY_DIR}/${OUTDIR}','${CMAKE_CURRENT_SOURCE_DIR}/${REFDIR}'))"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 60)  # Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endif()

endfunction()


function(compare_gemini_output TESTNAME TESTDIR REFDIR)

octave_compare(Compare2D ${TESTDIR} ${REFDIR})

matlab_compare(MatlabCompare2D ${TESTDIR} ${REFDIR})



endfunction()
