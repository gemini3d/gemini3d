
function(check_matlab_source_runs code)

if(NOT Matlab_MAIN_PROGRAM)
  set(ok false)
else()
  execute_process(COMMAND ${Matlab_MAIN_PROGRAM} -batch ${code}
    ERROR_QUIET OUTPUT_QUIET
    RESULT_VARIABLE ret
    TIMEOUT 60)  # Matlab takes a long time to start with lots of toolboxes
  if(ret EQUAL 0)
    set(MatlabOK true CACHE BOOL "Matlab is sufficiently new (>= R2019a) to run self-tests")
  else()
    set(MatlabOK false CACHE BOOL "Matlab is too old (< R2019a) to run self-tests")
  endif()
endif()

endfunction(check_matlab_source_runs)

find_package(Matlab COMPONENTS MAIN_PROGRAM)
check_matlab_source_runs("exit")