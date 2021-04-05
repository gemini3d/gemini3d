if(NOT matlab)
  return()
endif()

find_package(Matlab COMPONENTS MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  return()
endif()

cmake_minimum_required(VERSION 3.19...${CMAKE_VERSION})

FetchContent_Declare(MATGEMINI
GIT_REPOSITORY ${matgemini_git}
GIT_TAG ${matgemini_tag})

FetchContent_MakeAvailable(MATGEMINI)

if(NOT MATGEMINI_DIR)
  execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${matgemini_SOURCE_DIR}/setup.m'), gemini3d.fileio.expanduser('~');"
    RESULT_VARIABLE _ok
    TIMEOUT 90)

  if(_ok EQUAL 0)
    message(STATUS "MatGemini found: ${matgemini_SOURCE_DIR}")
    set(MATGEMINI_DIR ${matgemini_SOURCE_DIR} CACHE PATH "MatGemini path")
  endif()
endif()
