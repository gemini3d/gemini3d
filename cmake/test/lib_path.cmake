# Windows MSVC-based exe's (including Intel compiler on Windows) need DLL's on PATH.
if(CMAKE_VERSION VERSION_LESS 3.20)
  return()
endif()

if(NOT MSVC)
  message(DEBUG "SKIP: test_dll_path")
  return()
endif()

set(test_dll_path)

find_path(ZLIB_DLL NAMES zlib.dll
  NO_DEFAULT_PATH
  HINTS ${ZLIB_INCLUDE_DIR}/.. ${ZLIB_ROOT} ENV ZLIB_ROOT
  PATH_SUFFIXES bin
  DOC "MUMPS common header")

if(NOT ZLIB_DLL)
  set(ZLIB_DLL ${ZLIB_ROOT})
endif()

message(STATUS "${ZLIB_DLL}")

if(ZLIB_DLL)
  list(APPEND test_dll_path ${ZLIB_DLL})
endif()

cmake_path(APPEND_STRING test_dll_path ";$ENV{PATH}")
cmake_path(CONVERT "${test_dll_path}" TO_NATIVE_PATH_LIST test_dll_path NORMALIZE)

# this is the vital line, without it CMake set_tests_properties mangles the ENVIRONMENT
string(REPLACE ";" "\\;" test_dll_path "${test_dll_path}")

message(DEBUG "test_dll_path: ${test_dll_path}")
