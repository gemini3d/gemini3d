if(matlab AND NOT DEFINED MatlabOK)

find_package(Matlab COMPONENTS MAIN_PROGRAM)

# check if Matlab >= R2019 -batch is working
execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch exit
  OUTPUT_QUIET
  RESULT_VARIABLE ret
  TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes

if(ret EQUAL 0)
  set(MatlabOK true CACHE BOOL "Matlab can run self-tests")
else()
  set(MatlabOK false CACHE BOOL "Matlab is too old (< R2019a) to run self-tests")
endif()

endif()