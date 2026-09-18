cmake_minimum_required(VERSION 3.29)

if(NOT DEFINED file)
  message(FATAL_ERROR "Usage: cmake -Dfile=<file> -P ${CMAKE_CURRENT_LIST_FILE}")
endif()

if(IS_EXECUTABLE ${file})
  message(STATUS "${file} is executable")
else()
  message(STATUS "${file} is not executable")
endif()
