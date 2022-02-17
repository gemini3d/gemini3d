# Windows MSVC-based exe's (including Intel compiler on Windows) need DLL's on PATH.
# NOTE: CMake 3.22 added test property ENVIRONMENT_MODIFICATION that may do this more smoothly:
# https://cmake.org/cmake/help/git-stage/prop_test/ENVIRONMENT_MODIFICATION.html

if(NOT WIN32)
  return()
endif()

set(test_dll_path)

foreach(L ZLIB LAPACK)
  string(TOLOWER ${L} l)

  if(${L}_LIBRARIES)
    list(GET ${L}_LIBRARIES 0 t)
    cmake_path(GET t PARENT_PATH t)
  elseif(${L}_INCLUDE_DIR)
    list(GET ${L}_INCLUDE_DIR 0 t)
  else()
    continue()
  endif()

  cmake_path(GET t PARENT_PATH t2)

  find_path(${L}_DLL_DIR NAMES ${l}.dll lib${l}.dll
  NO_DEFAULT_PATH
  HINTS ${t2} ${${L}_ROOT} ENV CMAKE_PREFIX_PATH ENV ${L}_ROOT
  PATH_SUFFIXES bin
  DOC "${L} DLL PATH"
  )

  if(${L}_DLL_DIR)
    list(APPEND test_dll_path ${${L}_DLL_DIR})
  endif()
endforeach()

if(NOT test_dll_path)
  return()
endif()

cmake_path(APPEND_STRING test_dll_path ";$ENV{PATH}")
cmake_path(CONVERT "${test_dll_path}" TO_NATIVE_PATH_LIST test_dll_path NORMALIZE)

# this is the vital line, without it CMake set_tests_properties mangles the ENVIRONMENT
string(REPLACE ";" "\\;" test_dll_path "${test_dll_path}")

message(DEBUG "test_dll_path: ${test_dll_path}")
