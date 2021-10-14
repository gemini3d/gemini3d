# Windows MSVC-based exe's (including Intel compiler on Windows) need DLL's on PATH.
# NOTE: CMake 3.22 added test property ENVIRONMENT_MODIFICATION that may do this more smoothly:
# https://cmake.org/cmake/help/git-stage/prop_test/ENVIRONMENT_MODIFICATION.html

if(NOT MSVC)
  return()
endif()

set(test_dll_path)

find_path(ZLIB_DLL_DIR NAMES zlib.dll
  NO_DEFAULT_PATH
  HINTS ${ZLIB_INCLUDE_DIR}/.. ${ZLIB_ROOT} ENV ZLIB_ROOT
  PATH_SUFFIXES bin
  DOC "MUMPS common header")

if(NOT ZLIB_DLL_DIR)
  return()
endif()

set(test_dll_path ${ZLIB_DLL_DIR})

cmake_path(APPEND_STRING test_dll_path ";$ENV{PATH}")
cmake_path(CONVERT "${test_dll_path}" TO_NATIVE_PATH_LIST test_dll_path NORMALIZE)

# this is the vital line, without it CMake set_tests_properties mangles the ENVIRONMENT
string(REPLACE ";" "\\;" test_dll_path "${test_dll_path}")

message(DEBUG "test_dll_path: ${test_dll_path}")
