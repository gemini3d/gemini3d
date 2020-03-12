add_library(PDEparabolic OBJECT PDEparabolic.f90)
target_link_libraries(PDEparabolic PRIVATE const ${LAPACK_LIBRARIES})
target_include_directories(PDEparabolic PRIVATE ${LAPACK_INCLUDE_DIRS})
if(NOT LAPACK_LAPACK95_FOUND)
  target_sources(PDEparabolic PRIVATE ${PROJECT_SOURCE_DIR}/src/vendor/lapack95/gbsv.f90)
endif()

add_library(diffusion OBJECT diffusion.f90)
target_link_libraries(diffusion PRIVATE const grid mesh PDEparabolic)

if(NOT hdf5)
  return()
endif()

# --- tests

set(_tmpdifffn ${CMAKE_CURRENT_BINARY_DIR}/test_diffusion1d.h5)

add_executable(test_diffusion1d test_diffusion1D.f90)
# lapack needed for MacOS
target_link_libraries(test_diffusion1d PRIVATE PDEparabolic const
  ${LAPACK_LIBRARIES} h5fortran::h5fortran)

# WORKING_DIRECTORY is needed to write output file to repeatable path.
add_test(NAME unit:diffusion1 COMMAND $<TARGET_FILE:test_diffusion1d>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(unit:diffusion1 PROPERTIES
  TIMEOUT 5
  FIXTURES_SETUP GemDiff)
win32_env(unit:diffusion1)

add_test(NAME unit:diffusion1:python
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_diffusion1D.py ${_tmpdifffn})
set_tests_properties(unit:diffusion1:python PROPERTIES
  TIMEOUT 10
  REQUIRED_FILES ${_tmpdifffn}
  FIXTURES_REQUIRED GemDiff
  DISABLED ${python_disabled})
