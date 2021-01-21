if(NOT matlab)
  return()
endif()

find_package(Matlab COMPONENTS MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.19)
  message(FATAL_ERROR "MatGemini requires CMake >= 3.19")
endif()

FetchContent_Declare(matgemini
  GIT_REPOSITORY ${matgemini_git}
  GIT_TAG ${matgemini_tag})

FetchContent_MakeAvailable(matgemini)

if(NOT MATGEMINI_DIR)
  execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch "run('${CMAKE_CURRENT_SOURCE_DIR}/setup.m'), gemini3d.fileio.expanduser('~');"
    RESULT_VARIABLE _ok
    TIMEOUT 90)

  if(_ok EQUAL 0)
    message(STATUS "MatGemini found: ${matgemini_SOURCE_DIR}")
    set(MATGEMINI_DIR ${matgemini_SOURCE_DIR} CACHE PATH "MatGemini path")
  endif()
endif()
