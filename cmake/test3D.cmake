#---------------------------------
# --- testings
#----------------------------------
# Test data directories
set(TEST_REF_ROOT ${PROJECT_SOURCE_DIR}/../simulations)

# if the file version changes (new reference data uploaded) this URL changes as well
set(TEST_3D_URL https://zenodo.org/record/1473195/files/zenodo3d.zip?download=1)
set(TEST_3D_ARCHIVE ${TEST_REF_ROOT}/zenodo3d.zip)
set(TEST_3D_DIR ${TEST_REF_ROOT}/zenodo3d)
set(TEST_3D_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/test3d)


# --- ensure 3D reference data is available for self-test
if(NOT EXISTS ${TEST_3D_DIR})
  if(NOT EXISTS ${TEST_3D_ARCHIVE})
    file(DOWNLOAD ${TEST_3D_URL} ${TEST_3D_ARCHIVE} 
         SHOW_PROGRESS)
  endif()
  
  add_custom_command(TARGET gemini
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E tar -xf ${TEST_3D_ARCHIVE} 
                     WORKING_DIRECTORY ${TEST_REF_ROOT}
                     DEPENDS ${TEST_3D_ARCHIVE}
                     COMMENT "Unzipping 3D test files")
endif()


# --- test main exe

add_test(NAME Gemini3D 
         COMMAND mpirun -np 4 ./gemini ${CMAKE_SOURCE_DIR}/initialize/3Dtest/config.ini ${TEST_3D_OUTDIR} 
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(Gemini3D PROPERTIES 
                     TIMEOUT 900  # test should complete in under 15 minutes on ~ 2014 dual-core laptop
                     FIXTURES_SETUP GemMain)
                     

# --- evaluate output accuracy vs. reference from Matt's HPC
find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME Compare3D 
         COMMAND octave-cli -q --eval "compare_all('${TEST_3D_OUTDIR}','${TEST_3D_DIR}')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(Compare3D PROPERTIES
                     TIMEOUT 30)

endif()

# --- check that *.m syntax is Matlab compatible
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if (Matlab_FOUND)
add_test(NAME MatlabCompare3D
         COMMAND matlab -nojvm -r "exit(compare_all('${TEST_3D_OUTDIR}','${TEST_3D_DIR}'))"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
         
set_tests_properties(MatlabCompare3D PROPERTIES
                     TIMEOUT 60)  # Matlab with a lot of toolboxes takes ~ 15 seconds just to start

endif()
