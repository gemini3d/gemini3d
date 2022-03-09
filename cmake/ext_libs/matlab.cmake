include(FetchContent)

find_package(Matlab COMPONENTS MAIN_PROGRAM REQUIRED)
if(Matlab_VERSION_STRING AND Matlab_VERSION_STRING VERSION_LESS 9.9)
  message(STATUS "Matlab >= 9.9 required, found ${Matlab_VERSION_STRING}")
endif()

file(READ ${CMAKE_CURRENT_LIST_DIR}/../libraries.json _libj)
string(JSON matgemini_git GET ${_libj} matgemini git)

FetchContent_Declare(MATGEMINI
GIT_REPOSITORY ${matgemini_git}
GIT_TAG ${matgemini_tag}
INACTIVITY_TIMEOUT 15
)

FetchContent_Populate(MATGEMINI)

cmake_path(CONVERT "${matgemini_SOURCE_DIR};${matgemini_SOURCE_DIR}/matlab-stdlib/" TO_NATIVE_PATH_LIST MATLABPATH)

if(MATGEMINI_FOUND)
  return()
endif()

execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${matgemini_SOURCE_DIR}/setup.m'), stdlib.fileio.expanduser('~');"
TIMEOUT 90
RESULT_VARIABLE ret
ERROR_VARIABLE err
)

if(NOT ret EQUAL 0)
  message(FATAL_ERROR "MatGemini not available:
  ${ret}
  ${err}"
  )
endif()

set(MATGEMINI_FOUND true CACHE BOOL "MatGemini found")
