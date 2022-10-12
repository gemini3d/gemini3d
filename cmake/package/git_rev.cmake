# log the Git revision information for reproducibility
# This has the deficiency of not auto-updating on each build.
# CMake must be reconfigured to ensure the Git meta is updated.
# a possible workaround is in
# https://github.com/rpavlik/cmake-modules/blob/main/GetGitRevisionDescription.cmake

find_package(Git)

set(_max_len 80) # arbitrary limit, so as not to exceed maximum 132 character Fortran line length.
set(git_branch)
set(git_rev)
set(git_porcelain false)

if(NOT GIT_FOUND)
  return()
endif()

set(git_version ${GIT_VERSION_STRING})
string(SUBSTRING ${git_version} 0 ${_max_len} git_version)
# git branch --show-current requires Git >= 2.22, June 2019
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
OUTPUT_VARIABLE git_branch
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
if(_err EQUAL 0)
  string(SUBSTRING ${git_branch} 0 ${_max_len} git_branch)
else()
  set(git_branch)
endif()


# git describe --tags can make CI error fatal: No names found, cannot describe anything.
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
OUTPUT_VARIABLE git_rev
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
if(NOT _err EQUAL 0)
  # old Git
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE git_rev
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _err
  )
endif()
if(_err EQUAL 0)
  string(SUBSTRING ${git_rev} 0 ${_max_len} git_rev)
else()
  set(git_rev)
endif()
string(APPEND git_rev " ${PROJECT_VERSION}")

set(git_porcelain true)
execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
OUTPUT_VARIABLE _porcelain
OUTPUT_STRIP_TRAILING_WHITESPACE
RESULT_VARIABLE _err
)
if(_porcelain)
  set(git_porcelain false)
endif(_porcelain)
if(NOT _err EQUAL 0)
  set(git_porcelain false)
endif()

message(STATUS "${PROJECT_NAME}  git revision: ${git_rev}  git_branch: ${git_branch}  git_porcelain: ${git_porcelain}")
