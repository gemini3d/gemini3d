cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED exe)
  message(FATAL_ERROR "Please specify the executable to run -Dexe=<exe>")
endif()

if(UNIX)
  set(devnull INPUT_FILE /dev/null)
elseif(WIN32)
  set(devnull INPUT_FILE NUL)
endif()

execute_process(COMMAND ${exe}
${devnull}
COMMAND_ERROR_IS_FATAL ANY
)
