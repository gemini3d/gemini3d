function(clone_if_missing root_dir url)

if(ARGC GREATER 2)
  set(source_dir ${ARGV2})
else()
  set(source_dir ${root_dir})
endif()

if(IS_DIRECTORY ${source_dir})
  return()
endif()

find_package(Git REQUIRED)

execute_process(COMMAND ${GIT_EXECUTABLE} clone --depth 1 ${url} ${root_dir}
  RESULT_VARIABLE _gitstat
  TIMEOUT 60)
if(NOT _gitstat STREQUAL 0)
  message(FATAL_ERROR "could not Git clone ${url}, return code ${_gitstat}")
endif()

endfunction(clone_if_missing)
