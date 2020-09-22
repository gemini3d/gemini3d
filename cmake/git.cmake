function(git_download root_dir url tag)

find_package(Git REQUIRED)

execute_process(COMMAND ${GIT_EXECUTABLE} clone ${url} ${root_dir}
  RESULT_VARIABLE _gitstat
  TIMEOUT 120)
if(NOT _gitstat STREQUAL 0)
  message(FATAL_ERROR "could not Git clone ${url}, return code ${_gitstat}")
endif()

# use WORKING_DIRECTORY for legacy HPC Git
execute_process(COMMAND ${GIT_EXECUTABLE} checkout ${tag}
  RESULT_VARIABLE _gitstat
  TIMEOUT 30
  WORKING_DIRECTORY ${root_dir})
if(NOT _gitstat STREQUAL 0)
  message(FATAL_ERROR "could not Git checkout ${tag}, return code ${_gitstat}")
endif()

endfunction(git_download)
