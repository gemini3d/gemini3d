if(realbits EQUAL 32)
  add_library(mumps_interface mumps_real32.f90)
else()
  add_library(mumps_interface mumps_real64.f90)
endif()
target_include_directories(mumps_interface
  PUBLIC ${MUMPS_INCLUDE_DIRS}
  INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(mumps_interface PUBLIC ${MUMPS_LIBRARIES})
set_target_properties(mumps_interface PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(PDEelliptic PDEelliptic.f90)
target_link_libraries(PDEelliptic PRIVATE const mumps_interface MPI::MPI_Fortran)
target_include_directories(PDEelliptic INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(PDEelliptic PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

add_library(potential potential_comm_mumps.f90 potential_mumps.f90 boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PRIVATE calculus collision const grid mesh interp timeutils
  MPI::MPI_Fortran PDEelliptic mpimod reader)
target_include_directories(potential INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(potential PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

if(NOT hdf5)
  return()
endif()

set(_potential_testfile ${CMAKE_CURRENT_BINARY_DIR}/test_potential2d.h5)
add_executable(test_potential2d test_potential2D.f90)
target_link_libraries(test_potential2d const PDEelliptic MPI::MPI_Fortran h5fortran)
add_test(NAME potential2d COMMAND test_potential2d WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(potential2d PROPERTIES TIMEOUT 30
                     FIXTURES_SETUP GemPot)

if(OctaveOK)
add_test(NAME potential2d_octave
  COMMAND ${Octave_EXECUTABLE} -q --eval "test_potential2D('${_potential_testfile}')"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(potential2d_octave PROPERTIES
  REQUIRED_FILES ${_potential_testfile}
  TIMEOUT 15
  FIXTURES_REQUIRED GemPot)
endif(OctaveOK)

if(PythonOK)
add_test(NAME potential2d_python
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${_potential_testfile})
set_tests_properties(potential2d_python PROPERTIES
  REQUIRED_FILES ${_potential_testfile}
  TIMEOUT 10
  FIXTURES_REQUIRED GemPot)
endif(PythonOK)
