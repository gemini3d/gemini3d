include(${CMAKE_CURRENT_LIST_DIR}/git.cmake)

function(check_matgemini)

# other sources of info that aren't as easy to use
# * appdata\version.xml
# * appdata\prodcontents.json
# I have observed over the years that this directory name scheme is universally used.
string(REGEX MATCH "R([0-9][0-9][0-9][0-9][a-z])" Matlab_RELEASE ${Matlab_ROOT_DIR})
matlab_get_version_from_release_name(${Matlab_RELEASE} Matlab_VERSION)

if(Matlab_VERSION AND Matlab_VERSION VERSION_LESS 9.6)
  message(STATUS "Matlab >= R2019a required for -batch mode")
  return()
endif()

if(IS_DIRECTORY ${matgemini_dir})
  message(STATUS "Gemini-Matlab found: ${matgemini_dir}")
  set(matlab_ok true CACHE BOOL "Gemini-Matlab OK")
endif()

endfunction(check_matgemini)

# --- script
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  return()
endif()
# keep this in script so it's not scoped in function

set(matgemini_dir "${PROJECT_SOURCE_DIR}/../gemini-matlab/")
set(matgemini_url "https://github.com/gemini3d/gemini-matlab")

if(NOT matlab_ok)
  clone_if_missing(${matgemini_dir} ${matgemini_url})
  check_matgemini()
endif()
