add_library(PDEelliptic PDEelliptic.F90)
target_link_libraries(PDEelliptic PRIVATE const ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
target_compile_options(PDEelliptic PRIVATE ${FFLAGS})
target_include_directories(PDEelliptic PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_definitions(PDEelliptic PRIVATE REALBITS=${realbits})
set_target_properties(PDEelliptic PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_library(potential potential_comm_mumps.f90 potential_mumps.F90 boundary_conditions/potentialBCs_mumps.f90)
target_link_libraries(potential PRIVATE calculus collision const grid interp timeutils
  ${MUMPS_LIBRARIES} MPI::MPI_Fortran PDEelliptic mpimod)
target_compile_options(potential PRIVATE ${FFLAGS})
target_include_directories(potential PRIVATE ${MUMPS_INCLUDE_DIRS})
target_compile_definitions(potential PRIVATE REALBITS=${realbits})
set_target_properties(potential PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_potential2d test_potential2D.f90)
target_link_libraries(test_potential2d const PDEelliptic ${MUMPS_LIBRARIES} MPI::MPI_Fortran)
set_target_properties(test_potential2d PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME potential2d COMMAND test_potential2d WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(potential2d PROPERTIES TIMEOUT 30
                     FIXTURES_SETUP GemPot)

if(OctaveOK)
add_test(NAME OctavePotential2d
         COMMAND ${Octave_EXECUTABLE} -q --eval "test_potential2D('${CMAKE_BINARY_DIR}/test_potential2D.dat')"
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(OctavePotential2d  PROPERTIES
                     TIMEOUT 15
                     FIXTURES_REQUIRED GemPot)
endif(OctaveOK)

if(PythonOK)
add_test(NAME PythonPotential2d
  COMMAND ${Python3_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_potential2d.py ${CMAKE_BINARY_DIR}/test_potential2D.dat)
set_tests_properties(PythonPotential2d PROPERTIES
                     TIMEOUT 10
                     FIXTURES_REQUIRED GemPot)
endif(PythonOK)
