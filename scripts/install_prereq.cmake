# installs Gemini3D basic prereqs on Linux and MacOS, and Windows with MSYS2
# use by:
#
#  cmake -P scripts/install_prereq.cmake

cmake_minimum_required(VERSION 3.7...3.21)

if(WIN32)
  message(FATAL_ERROR "Please install Gemini prereqs on Windows via MSYS2 Terminal https://www.msys2.org/")
endif()

if(CMAKE_VERSION VERSION_LESS 3.20)
  message(FATAL_ERROR "Please first update CMake via
    cmake -P ${CMAKE_CURRENT_LIST_DIR}/install_cmake.cmake")
endif()

function(check_ninja)
  find_program(ninja NAMES ninja)
  if(ninja)
    execute_process(COMMAND ninja --version
      OUTPUT_STRIP_TRAILING_WHITESPACE
      TIMEOUT 5
      OUTPUT_VARIABLE ninja_ver)
    # don't fail output
  endif()

  if(NOT ninja OR ninja_ver VERSION_LESS 1.10.0)
    message(STATUS "Please install Ninja via
     cmake -P ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/install_ninja.cmake")
  endif()
endfunction(check_ninja)

# read prereqs

function(read_prereqs sys_id)

  file(READ ${CMAKE_CURRENT_LIST_DIR}/prereqs.json json)

  set(prereqs)
  string(JSON N LENGTH ${json} ${sys_id})
  math(EXPR N "${N}-1")
  foreach(i RANGE ${N})
    string(JSON _u GET ${json} ${sys_id} ${i})
    list(APPEND prereqs ${_u})
  endforeach()

  string(REPLACE ";" " " prereqs "${prereqs}")
  set(prereqs ${prereqs} PARENT_SCOPE)

endfunction(read_prereqs)

# detect platform

execute_process(COMMAND uname -s OUTPUT_VARIABLE id TIMEOUT 5)

if(id MATCHES "^MSYS")
  read_prereqs("msys2")
  execute_process(COMMAND pacman -S --needed ${prereqs})
  return()
endif()

if(APPLE)
  find_program(brew
    NAMES brew
    PATHS /usr/local /opt/homebrew
    PATH_SUFFIXES bin)

  if(NOT brew)
    message(FATAL_ERROR "We generally suggest installing Homebrew package manager https://brew.sh")
  endif()

  read_prereqs("brew")
  execute_process(COMMAND ${brew} install ${prereqs})
  return()
endif()

find_program(apt NAMES apt)
if(apt)
  read_prereqs("apt")
  execute_process(COMMAND apt install --no-install-recommends ${prereqs})
  # don't check return code as can be non-zero if packages already installed
  check_ninja()
  return()
endif()

find_program(yum NAMES yum)
if(yum)
  read_prereqs("yum")
  execute_process(COMMAND yum install ${prereqs})
  # don't check return code as can be non-zero if packages already installed
  check_ninja()
  return()
endif()

find_program(pacman NAMES pacman)
if(pacman)
  read_prereqs("pacman")
  execute_process(COMMAND pacman -S --needed ${prereqs})

  check_ninja()
  return()
endif()
