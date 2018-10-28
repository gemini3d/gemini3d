add_executable(test_diffusion1D diffusion/test_diffusion1D.f90)
target_link_libraries(test_diffusion1D diffusion const)

add_test(NAME diffusion1D COMMAND test_diffusion1D)
set_tests_properties(diffusion1D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemDiff)

find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)

add_test(NAME OctaveDiffusion1D 
         COMMAND octave-cli -q --eval "test_diffusion1D('${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffusion)

set_tests_properties(OctaveDiffusion1D  PROPERTIES
                     TIMEOUT 5
                     FIXTURES_REQUIRED GemDiff)
          
endif()
