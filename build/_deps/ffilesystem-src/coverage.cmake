# run like:
#
#  ctest -S coverage.cmake

list(APPEND opts -DCMAKE_BUILD_TYPE=Debug)

set(CTEST_TEST_TIMEOUT 60)
# takes effect only if test property TIMEOUT is not set

find_program(exe NAMES gcov REQUIRED)
set(CTEST_COVERAGE_COMMAND ${exe})

set(check_flags --coverage)

if(check_flags)
  list(APPEND opts
  -DCMAKE_C_FLAGS_DEBUG=${check_flags}
  -DCMAKE_CXX_FLAGS_DEBUG=${check_flags}
  -DCMAKE_EXE_LINKER_FLAGS_INIT=${check_flags}
  )
endif()

set(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set(CTEST_BINARY_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/build-${CTEST_MEMORYCHECK_TYPE})
set(CTEST_BUILD_CONFIGURATION Debug)

if(DEFINED ENV{CMAKE_GENERATOR})
  set(CTEST_CMAKE_GENERATOR $ENV{CMAKE_GENERATOR})
else()
  set(CTEST_CMAKE_GENERATOR "Ninja")
endif()
set(CTEST_BUILD_FLAGS -j)

message(STATUS "Checker ${CTEST_MEMORYCHECK_TYPE}: ${CTEST_MEMORYCHECK_COMMAND}")

ctest_start(Experimental)

ctest_configure(
OPTIONS "${opts}"
RETURN_VALUE ret
CAPTURE_CMAKE_ERROR err
)
if(NOT (ret EQUAL 0 AND err EQUAL 0))
  message(FATAL_ERROR "CMake configure failed:  ${ret}   ${err}")
endif()

ctest_build(
RETURN_VALUE ret
CAPTURE_CMAKE_ERROR err
)
if(NOT (ret EQUAL 0 AND err EQUAL 0))
  message(FATAL_ERROR "CMake build failed:  ${ret}   ${err}")
endif()

ctest_coverage(
RETURN_VALUE ret
CAPTURE_CMAKE_ERROR err
)

if(NOT (ret EQUAL 0 AND err EQUAL 0))
  message(FATAL_ERROR "Coverage check failed:  ${ret}   ${err}")
endif()
