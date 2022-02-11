include(FetchContent)

find_package(Matlab 9.9 COMPONENTS MAIN_PROGRAM REQUIRED)

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
