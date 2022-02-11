# prints Gemini3D prereqs on stdout
#  cmake -P scripts/requirements.cmake

cmake_minimum_required(VERSION 3.20...3.23)

cmake_path(SET prereq_file NORMALIZE ${CMAKE_CURRENT_LIST_DIR}/../requirements.json)

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

execute_process(COMMAND uname -s OUTPUT_VARIABLE id TIMEOUT 5)

if(id MATCHES "^MINGW64")
  read_prereqs("msys2")
  execute_process(COMMAND ${CMAKE_COMMAND} -E echo "${cmd} ${prereqs}" TIMEOUT 2)
  return()
endif()

if(APPLE)
  set(names brew port)
elseif(UNIX)
  set(names apt yum pacman zypper)
elseif(WIN32)
  message(FATAL_ERROR "Windows: suggest MSYS2 or WSL")
endif()

foreach(t IN LISTS names)
  find_program(${t} NAMES ${t})
  if(${t})
    read_prereqs(${t})
    execute_process(COMMAND ${CMAKE_COMMAND} -E echo "${cmd} ${prereqs}" TIMEOUT 2)
    return()
  endif()
endforeach()

message(FATAL_ERROR "Package manager not found ${names}")
