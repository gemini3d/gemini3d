# This CTest script assumes you've already configured and built, and just want to submit a CDash result
# run by:
# ctest -S cdash.cmake

set(CTEST_BUILD_CONFIGURATION "Release")
set(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
if(NOT DEFINED CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/build)
endif()

set(CTEST_CMAKE_GENERATOR "Ninja")

site_name(CTEST_SITE)
set(CTEST_BUILD_NAME ${CMAKE_SYSTEM_NAME})

ctest_start("Experimental" ${CTEST_SOURCE_DIRECTORY} ${CTEST_BINARY_DIRECTORY})
# ctest_configure(BUILD ${CTEST_BINARY_DIRECTORY} SOURCE ${CTEST_SOURCE_DIRECTORY})
# ctest_build(BUILD ${CTEST_BINARY_DIRECTORY} CONFIGURATION ${CTEST_BUILD_CONFIGURATION})
ctest_test(BUILD ${CTEST_BINARY_DIRECTORY})
ctest_submit()