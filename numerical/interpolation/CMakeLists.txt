add_library(interp interpolation.f90)
target_link_libraries(interp PRIVATE const)
target_compile_options(interp PRIVATE ${FFLAGS})
set_target_properties(interp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_interp1 testinterp1.f90)
target_link_libraries(test_interp1 PRIVATE const interp)
target_compile_options(test_interp1 PRIVATE ${FFLAGS})
set_target_properties(test_interp1 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_interp2 testinterp2.f90)
target_link_libraries(test_interp2 PRIVATE const interp)
target_compile_options(test_interp2 PRIVATE ${FFLAGS})
set_target_properties(test_interp2 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_interp3 testinterp3.f90)
target_link_libraries(test_interp3 PRIVATE const interp)
target_compile_options(test_interp3 PRIVATE ${FFLAGS})
set_target_properties(test_interp3 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})


# WORKING_DIRECTORY is needed to write output file to repeatable path.
add_test(NAME Interp1d COMMAND test_interp1 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(Interp1d PROPERTIES
                     TIMEOUT 10
                     FIXTURES_SETUP GemInterp)

add_test(NAME Interp2d COMMAND test_interp2 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(Interp2d PROPERTIES
                     TIMEOUT 10
                     FIXTURES_SETUP GemInterp)

add_test(NAME Interp3d COMMAND test_interp3 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(Interp3d PROPERTIES
                     TIMEOUT 30  # need 30 for slow CPU
                     FIXTURES_SETUP GemInterp)

if(PythonOK)
add_test(NAME InterpPython
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testinterp.py ${CMAKE_BINARY_DIR}/output2D.dat)
set_tests_properties(InterpPython PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemInterp)

add_test(NAME Interp3Python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testinterp3.py ${CMAKE_BINARY_DIR}/output3D.dat $<TARGET_FILE:test_interp3>)
set_tests_properties(Interp3Python PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemInterp)
endif(PythonOK)

if(OctaveOK)
add_test(NAME InterpOctave
         COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp('${CMAKE_BINARY_DIR}/output2D.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(InterpOctave PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemInterp)

add_test(NAME Interp3Octave
         COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp3('${CMAKE_BINARY_DIR}/output3D.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(Interp3Octave PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemInterp)
endif(OctaveOK)
