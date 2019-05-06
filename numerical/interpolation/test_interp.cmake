add_executable(test_interp1 interpolation/testinterp1.f90)
target_link_libraries(test_interp1 interp)
target_compile_options(test_interp1 PRIVATE ${FFLAGS})
set_target_properties(test_interp1 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})


add_executable(test_interp2 interpolation/testinterp2.f90)
target_link_libraries(test_interp2 interp)
target_compile_options(test_interp2 PRIVATE ${FFLAGS})
set_target_properties(test_interp2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})


add_test(NAME Interp1D COMMAND test_interp1)
set_tests_properties(Interp1D PROPERTIES
                     TIMEOUT 5)

add_test(NAME Interp2D COMMAND test_interp2)
set_tests_properties(Interp2D PROPERTIES
                     TIMEOUT 5
                     FIXTURES_SETUP GemInterp)


if(OctaveOK)

add_test(NAME OctaveInterp
         COMMAND ${Octave_EXECUTABLE} --eval "testinterp('${CMAKE_CURRENT_BINARY_DIR}/output2D.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


set_tests_properties(OctaveInterp PROPERTIES
                     TIMEOUT 5
                     FIXTURES_REQUIRED GemInterp)

endif()
