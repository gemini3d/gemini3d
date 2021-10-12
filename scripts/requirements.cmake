# prints Gemini3D prereqs on stderr
#  cmake -P scripts/requirements.cmake

cmake_minimum_required(VERSION 3.7...3.22)

set(prereq_file ${CMAKE_CURRENT_LIST_DIR}/../requirements.json)
cmake_path(ABSOLUTE_PATH prereq_file NORMALIZE)

# --- helper functions

function(read_prereqs sys_id)

  file(READ ${prereq_file} json)

  set(prereqs)
  string(JSON N LENGTH ${json} ${sys_id} pkgs)
  math(EXPR N "${N}-1")
  foreach(i RANGE ${N})
    string(JSON _u GET ${json} ${sys_id} pkgs ${i})
    list(APPEND prereqs ${_u})
  endforeach()
  string(REPLACE ";" " " prereqs "${prereqs}")

  string(JSON cmd GET ${json} ${sys_id} cmd)

  set(prereqs ${prereqs} PARENT_SCOPE)
  set(cmd ${cmd} PARENT_SCOPE)

endfunction(read_prereqs)

# --- main program

if(CMAKE_VERSION VERSION_LESS 3.20)
  message(NOTICE "cmake -P ${CMAKE_CURRENT_LIST_DIR}/install_cmake.cmake")
endif()

execute_process(COMMAND uname -s OUTPUT_VARIABLE id TIMEOUT 5)

if(APPLE)
  # avoids MacOS-only tools with conflicting names
  find_program(brew NAMES brew)
  find_program(port NAMES port)
else()
  find_program(apt NAMES apt)
  find_program(yum NAMES yum)
  find_program(pacman NAMES pacman)
endif()

if(apt)
  read_prereqs("apt")
elseif(yum)
  read_prereqs("yum")
elseif(pacman)
  read_prereqs("pacman")
elseif(brew)
  read_prereqs("brew")
elseif(port)
  read_prereqs("port")
elseif(id MATCHES "^MSYS")
  read_prereqs("msys2")
endif()

message(NOTICE "${cmd} ${prereqs}")
