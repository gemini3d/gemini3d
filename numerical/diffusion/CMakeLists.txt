add_library(PDEparabolic PDEparabolic.f90)
target_link_libraries(PDEparabolic PRIVATE const ${LAPACK_LIBRARIES})
target_include_directories(PDEparabolic PRIVATE ${LAPACK_INCLUDE_DIRS})
if(NOT LAPACK_LAPACK95_FOUND)
  target_sources(PDEparabolic PRIVATE ${CMAKE_SOURCE_DIR}/vendor/lapack95/gbsv.f90)
endif()

add_library(diffusion diffusion.f90)
target_link_libraries(diffusion PRIVATE const grid mesh PDEparabolic)


add_executable(test_diffusion1D test_diffusion1D.f90)
# lapack for mac
target_link_libraries(test_diffusion1D PRIVATE PDEparabolic const ${LAPACK_LIBRARIES})

# WORKING_DIRECTORY is needed to write output file to repeatable path.
add_test(NAME diffusion1D COMMAND test_diffusion1D WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(diffusion1D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemDiff)


if(OctaveOK)
add_test(NAME OctaveDiffusion1D
         COMMAND ${Octave_EXECUTABLE} -q --eval "test_diffusion1D('${CMAKE_BINARY_DIR}/test_diffusion1d.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(OctaveDiffusion1D  PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemDiff)
endif(OctaveOK)

if(PythonOK)
add_test(NAME PythonDiffusion1D
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion1D.py ${CMAKE_BINARY_DIR}/test_diffusion1d.dat)
set_tests_properties(PythonDiffusion1D PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemDiff)
endif(PythonOK)
