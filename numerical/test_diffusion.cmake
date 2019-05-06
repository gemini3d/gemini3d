add_executable(test_diffusion1D diffusion/test_diffusion1D.f90)
target_link_libraries(test_diffusion1D diffusion const)
set_target_properties(test_diffusion1D PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME diffusion1D COMMAND test_diffusion1D)
set_tests_properties(diffusion1D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemDiff)


if(OctaveOK)

add_test(NAME OctaveDiffusion1D
         COMMAND ${Octave_EXECUTABLE} --eval "test_diffusion1D('${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/diffusion)

set_tests_properties(OctaveDiffusion1D  PROPERTIES
                     TIMEOUT 5
                     FIXTURES_REQUIRED GemDiff)

endif()
