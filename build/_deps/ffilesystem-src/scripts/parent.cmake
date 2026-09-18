# compare path behavior between
#
# CMake: cmake_path() and get_filename_component()
# Python: pathlib.Path()
# Ffilesystem
#

cmake_minimum_required(VERSION 3.20)

find_package(Python COMPONENTS Interpreter)

find_program(fscli NAMES filesystem_cli PATHS ${CMAKE_CURRENT_LIST_DIR}/../build HINTS ${ffilesystem_ROOT})
if(fscli)
 execute_process(COMMAND ${fscli} cpp OUTPUT_VARIABLE ffs_cpp OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

# find_package(Matlab COMPONENTS MAIN_PROGRAM)  # requires Project role due to add_library()
find_program(Matlab_EXECUTABLE NAMES matlab HINTS $ENV{Matlab_ROOT})

# set(CMAKE_EXECUTE_PROCESS_COMMAND_ECHO STDOUT)



function(normal)

foreach(in IN ITEMS a/b/ a/b/.)

# cmake_path(NORMAL_PATH in)

cmake_path(GET in PARENT_PATH out)
message(STATUS "CMake: cmake_path(GET ${in} PARENT_PATH) => ${out}")

get_filename_component(out ${in} DIRECTORY)
message(STATUS "CMake: get_filename_component(${in} DIRECTORY) => ${out}")

if(Python_Interpreter_FOUND)

execute_process(
  COMMAND ${Python_EXECUTABLE} -c "from pathlib import Path; print(Path('${in}').parent)"
  OUTPUT_VARIABLE py OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE ret
)
if(ret EQUAL 0)
message(STATUS "Python ${Python_VERSION}: pathlib.Path('${in}').parent => ${py}")
endif()

endif(Python_Interpreter_FOUND)

if(fscli)

execute_process(
  COMMAND ${fscli} "parent" "${in}"
  OUTPUT_VARIABLE ffs OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE ret
)
if(ret EQUAL 0)
message(STATUS "Ffilesystem C++ ${ffs_cpp}: fs_parent(${in}) => ${ffs}")
endif()

endif(fscli)

message("")


endforeach()

endfunction(normal)


normal()
