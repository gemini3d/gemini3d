if(matlab)
  if(NOT DEFINED matlab_disabled)
    find_package(Matlab COMPONENTS MAIN_PROGRAM)

    # check if Matlab >= R2019 -batch is working
    execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch exit
      OUTPUT_QUIET
      RESULT_VARIABLE _ok
      TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes

    if(_ok EQUAL 0)
      set(matlab_disabled false CACHE BOOL "0: Matlab new enough for -batch")
    else()
      set(matlab_disabled true)
    endif()
  endif()
else()
  set(matlab_disabled true)
endif(matlab)