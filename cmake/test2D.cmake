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
         SHOW_PROGRESS)
  endif()
  
  add_custom_command(TARGET gemini
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E tar -xf ${TEST_2D_ARCHIVE} 
                     WORKING_DIRECTORY ${TEST_REF_ROOT}
                     DEPENDS ${TEST_2D_ARCHIVE}
                     COMMENT "Unzipping 2D test files")
endif()


# --- test main exe

add_test(NAME Gemini2D 
         COMMAND mpirun -np 4 ./gemini ${CMAKE_SOURCE_DIR}/initialize/2Dtest/config.ini ${TEST_2D_OUTDIR} 
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(Gemini2D PROPERTIES 
                     TIMEOUT 900  # test should complete in under 15 minutes on ~ 2014 dual-core laptop
                     FIXTURES_SETUP GemMain)
                     

# --- evaluate output accuracy vs. reference from Matt's HPC
find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME ReferenceCompare 
         COMMAND octave-cli -q --eval "compare_all('${TEST_2D_OUTDIR}','${TEST_2D_DIR}')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(ReferenceCompare PROPERTIES
                     TIMEOUT 5)

endif()

# --- check that *.m syntax is Matlab compatible
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if (Matlab_FOUND)
add_test(NAME MatlabCompare 
         COMMAND matlab -nojvm -r "exit(compare_all('${TEST_2D_OUTDIR}','${TEST_2D_DIR}'))"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(MatlabCompare PROPERTIES
                     TIMEOUT 30)  # Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endif()
