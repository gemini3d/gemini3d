#---------------------------------
# --- testings
#----------------------------------
# Test data directories
set(TEST_REF_ROOT ${PROJECT_SOURCE_DIR}/../simulations)

# if the file version changes (new reference data uploaded) this URL changes as well
set(TEST_2D_URL https://zenodo.org/record/1470073/files/zenodo2D.zip?download=1)
set(TEST_2D_ARCHIVE ${TEST_REF_ROOT}/zenodo2D.zip)
set(TEST_2D_DIR ${TEST_REF_ROOT}/zenodo2D)
set(TEST_2D_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/test2d)


# --- ensure 2D reference data is available for self-test
if(NOT EXISTS ${TEST_2D_DIR})
  if(NOT EXISTS ${TEST_2D_ARCHIVE})
    file(DOWNLOAD ${TEST_2D_URL} ${TEST_2D_ARCHIVE} 
         SHOW_PROGRESS
		 EXPECTED_HASH MD5=7463818a6e96ae8f5fff07dba2e3cdb8)
  endif()
  
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf ${TEST_2D_ARCHIVE} 
                  WORKING_DIRECTORY ${TEST_REF_ROOT})
endif()


# --- test main exe
cmake_host_system_information(RESULT PHYSRAM QUERY AVAILABLE_PHYSICAL_MEMORY)
if(${PHYSRAM} LESS 1000)
  set(NTEST 1)
else()
  set(NTEST 4)
endif()

add_test(NAME Gemini2D 
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NTEST} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gemini ${CMAKE_SOURCE_DIR}/initialize/2Dtest/config.ini ${TEST_2D_OUTDIR} 
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(Gemini2D PROPERTIES 
                     TIMEOUT 900  # test should complete in under 15 minutes on ~ 2014 dual-core laptop
                     FIXTURES_SETUP GemMain)
                     

# --- evaluate output accuracy vs. reference from Matt's HPC
find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME Compare2D 
         COMMAND ${OCTAVE_EXECUTABLE} -q --eval "compare_all('${TEST_2D_OUTDIR}','${TEST_2D_DIR}')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(Compare2D PROPERTIES
                     TIMEOUT 30)

endif()

# --- check that *.m syntax is Matlab compatible
find_package(Matlab QUIET COMPONENTS MAIN_PROGRAM)
if (Matlab_FOUND)
add_test(NAME MatlabCompare2D
         COMMAND ${Matlab_MAIN_PROGRAM} -nojvm -r "exit(compare_all('${TEST_2D_OUTDIR}','${TEST_2D_DIR}'))"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(MatlabCompare2D PROPERTIES
                     TIMEOUT 60)  # Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endif()
