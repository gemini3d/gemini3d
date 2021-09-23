# --- detect Homebrew / MacPorts and hint location
# this helps avoid issues with Anaconda overriding HDF5 with its broken compiler wrapper
if(NOT APPLE)
  return()
endif()

if(NOT DEFINED ENV{HOMEBREW_PREFIX})
  find_program(HOMEBREW NAMES brew)
  if(HOMEBREW)
    cmake_path(GET HOMEBREW PARENT_PATH HOMEBREW_PREFIX)
  endif()
endif()

if(NOT DEFINED ENV{MACPORTS_PREFIX})
  find_program(MACPORTS NAMES port)
  if(MACPORTS)
    cmake_path(GET MACPORTS PARENT_PATH MACPORTS_PREFIX)
  endif()
endif()
