include(${CMAKE_CURRENT_LIST_DIR}/git.cmake)

function(check_matgemini)

# other sources of info that aren't as easy to use
# * appdata\version.xml
# * appdata\prodcontents.json
# I have observed over the years that this directory name scheme is universally used.
# NOTE: some do not put an "R" in front that is /opt/matlab/2020a/... instead of /opt/matlab/R2020a/...
#
set(Matlab_VERSION)
string(REGEX MATCH "([0-9][0-9][0-9][0-9][a-z])" Matlab_RELEASE ${Matlab_ROOT_DIR})
if(Matlab_RELEASE) # don't fail if installed in custom directory
  matlab_get_version_from_release_name(${Matlab_RELEASE} Matlab_VERSION)
endif()

if(Matlab_VERSION AND Matlab_VERSION VERSION_LESS 9.6)
  message(STATUS "Matlab >= R2019a required for -batch mode")
  return()
endif()

if(IS_DIRECTORY ${matgemini_dir})
  message(STATUS "MatGemini found: ${matgemini_dir}")
  set(matlab_ok true CACHE BOOL "MatGemini OK")
endif()

endfunction(check_matgemini)

# --- script
find_package(Matlab COMPONENTS MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  return()
endif()
# keep this in script so it's not scoped in function

set(matgemini_dir "${PROJECT_SOURCE_DIR}/../mat_gemini/")
set(matgemini_url "https://github.com/gemini3d/mat_gemini")

if(NOT matlab_ok)
  clone_if_missing(${matgemini_dir} ${matgemini_url})
  check_matgemini()
endif()
