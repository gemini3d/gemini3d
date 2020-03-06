if(matlab)

if(DEFINED matlab_disabled)
  return()
endif()

find_package(Matlab COMPONENTS MAIN_PROGRAM)

# check if Matlab >= R2019 -batch is working
execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch exit
  OUTPUT_QUIET
  RESULT_VARIABLE matlab_disabled
  TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes

set(matlab_disabled ${matlab_disabled} CACHE STRING "0: Matlab new enough for -batch")
else()
set(matlab_disabled true CACHE BOOL "Matlab new enough for -batch")
endif()