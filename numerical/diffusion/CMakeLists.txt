add_library(diffusion diffusion.F90)
target_link_libraries(diffusion PRIVATE const grid ${LAPACK_LIBRARIES})
target_compile_options(diffusion PRIVATE ${FFLAGS})
target_include_directories(diffusion PRIVATE ${LAPACK_INCLUDE_DIRS})
if(NOT LAPACK_LAPACK95_FOUND)
  target_sources(diffusion PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/lapack95/gbsv.F90)
  target_compile_definitions(diffusion PRIVATE REALBITS=${realbits})
endif()
set_target_properties(diffusion PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})


add_executable(test_diffusion1D test_diffusion1D.f90)
target_link_libraries(test_diffusion1D diffusion const)
set_target_properties(test_diffusion1D PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME diffusion1D COMMAND test_diffusion1D)
set_tests_properties(diffusion1D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemDiff)


if(OctaveOK)
add_test(NAME OctaveDiffusion1D
         COMMAND ${Octave_EXECUTABLE} --eval "test_diffusion1D('${CMAKE_BINARY_DIR}/test_diffusion1d.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(OctaveDiffusion1D  PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemDiff)
endif(OctaveOK)
