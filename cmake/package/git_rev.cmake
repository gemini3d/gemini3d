# log the Git revision information for reproducibility
# This has the deficiency of not auto-updating on each build.
# CMake must be reconfigured to ensure the Git meta is updated.
# a possible workaround is in
# https://github.com/rpavlik/cmake-modules/blob/main/GetGitRevisionDescription.cmake

find_package(Git)

set(_max_len 80) # arbitrary limit, so as not to exceed maximum 132 character Fortran line length.
set(git_branch)
set(git_rev)
set(git_porcelain)
set(git_remote)

if(NOT GIT_FOUND)
  return()
endif()

set(git_version ${GIT_VERSION_STRING})
string(SUBSTRING ${git_version} 0 ${_max_len} git_version)
# git branch --show-current requires Git >= 2.22, June 2019
execute_process(COMMAND ${GIT_EXECUTABLE} -C ${PROJECT_SOURCE_DIR} rev-parse --abbrev-ref HEAD
OUTPUT_VARIABLE git_branch
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
if(_err EQUAL 0)
  string(SUBSTRING ${git_branch} 0 ${_max_len} git_branch)
else()
  set(git_branch)
endif()


execute_process(COMMAND ${GIT_EXECUTABLE} -C ${PROJECT_SOURCE_DIR} describe --tags
OUTPUT_VARIABLE git_rev
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
if(_err EQUAL 0)
  string(SUBSTRING ${git_rev} 0 ${_max_len} git_rev)
else()
  set(git_rev)
endif()
string(APPEND git_rev " ${PROJECT_VERSION}")

execute_process(COMMAND ${GIT_EXECUTABLE} -C ${PROJECT_SOURCE_DIR} status --porcelain
OUTPUT_VARIABLE _porcelain
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
string(LENGTH "${_porcelain}" _L)
if(_L EQUAL 0 AND _err EQUAL 0)
  set(git_porcelain "clean")
else()
  set(git_porcelain "dirty")
endif()

set(git_origin "origin")
# Git default remote name

execute_process(COMMAND ${GIT_EXECUTABLE} -C ${PROJECT_SOURCE_DIR} remote get-url ${git_origin}
OUTPUT_VARIABLE git_remote
OUTPUT_STRIP_TRAILING_WHITESPACE
)
# allow error message on stdout to propagate to CMake output

message(STATUS "${PROJECT_NAME} ${git_remote} ${git_rev} ${git_branch} ${git_porcelain}")
