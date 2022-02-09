include(FetchContent)

find_package(Matlab COMPONENTS MAIN_PROGRAM REQUIRED)

FetchContent_Declare(MATGEMINI
GIT_REPOSITORY ${matgemini_git}
GIT_TAG ${matgemini_tag}
INACTIVITY_TIMEOUT 15
)

FetchContent_MakeAvailable(MATGEMINI)

cmake_path(CONVERT "${matgemini_SOURCE_DIR};${matgemini_SOURCE_DIR}/matlab-stdlib/" TO_NATIVE_PATH_LIST MATLABPATH)

if(MATGEMINI_DIR)
  return()
endif()

execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${matgemini_SOURCE_DIR}/setup.m'), stdlib.fileio.expanduser('~');"
TIMEOUT 90
COMMAND_ERROR_IS_FATAL ANY
)

message(STATUS "MatGemini found: ${matgemini_SOURCE_DIR}")
set(MATGEMINI_DIR ${matgemini_SOURCE_DIR} CACHE PATH "MatGemini path")
