cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED exe)
  message(FATAL_ERROR "Please specify the executable to run -Dexe=<exe>")
endif()

execute_process(COMMAND ${exe}
INPUT_FILE ${CMAKE_CURRENT_LIST_DIR}/cli_exercise.txt
COMMAND_ERROR_IS_FATAL ANY
COMMAND_ECHO STDOUT
OUTPUT_VARIABLE out
ERROR_VARIABLE err
)

if(out MATCHES "TRACE")
  message(FATAL_ERROR "TRACE found in stdout")
endif()

if(err MATCHES "TRACE")
  message(FATAL_ERROR "TRACE found in stderr")
endif()

if(NOT err STREQUAL "")
  message(FATAL_ERROR "stderr output is not empty:
  ${err}")
endif()

message(STATUS "stdout:
${out}")
