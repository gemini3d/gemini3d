add_library(interp interpolation.f90)
target_link_libraries(interp PRIVATE const)

# --- tests
add_executable(test_interp1 testinterp1.f90)
target_link_libraries(test_interp1 PRIVATE const interp)
add_test(NAME Interp1d COMMAND test_interp1 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(Interp1d PROPERTIES
  TIMEOUT 10
  FIXTURES_SETUP GemInterp)

if(NOT hdf5)
  return()
endif()

set(_testfile_interp2 ${CMAKE_CURRENT_BINARY_DIR}/output2d.h5)
set(_testfile_interp3 ${CMAKE_CURRENT_BINARY_DIR}/output3d.h5)

add_executable(test_interp2 testinterp2.f90)
target_link_libraries(test_interp2 PRIVATE interp const h5fortran)
add_test(NAME Interp2d COMMAND test_interp2 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(Interp2d PROPERTIES
  TIMEOUT 10
  FIXTURES_SETUP GemInterp)

add_executable(test_interp3 testinterp3.f90)
target_link_libraries(test_interp3 PRIVATE const interp h5fortran)
add_test(NAME Interp3d COMMAND test_interp3 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(Interp3d PROPERTIES
  TIMEOUT 30  # need 30 for slow CPU
  FIXTURES_SETUP GemInterp)


if(PythonOK)
add_test(NAME Interp2Python
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testinterp.py ${_testfile_interp2})
set_tests_properties(Interp2Python PROPERTIES
  REQUIRED_FILES ${_testfile_interp2}
  TIMEOUT 10
  FIXTURES_REQUIRED GemInterp)

add_test(NAME Interp3Python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testinterp3.py ${_testfile_interp3})
set_tests_properties(Interp3Python PROPERTIES
  REQUIRED_FILES ${_testfile_interp3}
  TIMEOUT 10
  FIXTURES_REQUIRED GemInterp)
endif(PythonOK)

if(OctaveOK)
add_test(NAME Interp2Octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp('${_testfile_interp2}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(Interp2Octave PROPERTIES
  REQUIRED_FILES ${_testfile_interp2}
  TIMEOUT 10
  FIXTURES_REQUIRED GemInterp)

add_test(NAME Interp3Octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp3('${_testfile_interp3}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(Interp3Octave PROPERTIES
  REQUIRED_FILES ${_testfile_interp3}
  TIMEOUT 10
  FIXTURES_REQUIRED GemInterp)
endif(OctaveOK)
