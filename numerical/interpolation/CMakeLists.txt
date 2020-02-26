add_library(interp OBJECT interpolation.f90)
target_link_libraries(interp PRIVATE const)

# --- tests
add_executable(test_interp1 testinterp1.f90)
target_link_libraries(test_interp1 PRIVATE const interp)
add_test(NAME unit:interp1d COMMAND test_interp1 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:interp1d PROPERTIES
  TIMEOUT 10
  FIXTURES_SETUP GemInterp)

if(NOT hdf5)
  return()
endif()

set(_testfile_interp2 ${CMAKE_CURRENT_BINARY_DIR}/output2d.h5)
set(_testfile_interp3 ${CMAKE_CURRENT_BINARY_DIR}/output3d.h5)

add_executable(test_interp2 testinterp2.f90)
target_link_libraries(test_interp2 PRIVATE interp const h5fortran::h5fortran)
add_test(NAME unit:interp2d COMMAND test_interp2 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:interp2d PROPERTIES
  TIMEOUT 15
  FIXTURES_SETUP GemInterp)

add_executable(test_interp3 testinterp3.f90)
target_link_libraries(test_interp3 PRIVATE const interp h5fortran::h5fortran)
add_test(NAME unit:interp3d COMMAND test_interp3 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:interp3d PROPERTIES
  TIMEOUT 30  # need 30 for slow CPU
  FIXTURES_SETUP GemInterp)


if(PythonOK)
add_test(NAME unit:interp2:python
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testinterp.py ${_testfile_interp2})
set_tests_properties(unit:interp2:python PROPERTIES
  REQUIRED_FILES ${_testfile_interp2}
  TIMEOUT 30
  FIXTURES_REQUIRED GemInterp)

add_test(NAME unit:interp3:python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testinterp3.py ${_testfile_interp3})
set_tests_properties(unit:interp3:python PROPERTIES
  REQUIRED_FILES ${_testfile_interp3}
  TIMEOUT 30
  RUN_SERIAL true  # due to memory use
  FIXTURES_REQUIRED GemInterp)
endif(PythonOK)

if(OctaveOK)
add_test(NAME unit:interp2:octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp('${_testfile_interp2}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(unit:interp2:octave PROPERTIES
  REQUIRED_FILES ${_testfile_interp2}
  TIMEOUT 30
  FIXTURES_REQUIRED GemInterp)

add_test(NAME unit:interp3:octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "testinterp3('${_testfile_interp3}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(unit:interp3:octave PROPERTIES
  REQUIRED_FILES ${_testfile_interp3}
  RUN_SERIAL true
  TIMEOUT 30
  FIXTURES_REQUIRED GemInterp)
endif(OctaveOK)
