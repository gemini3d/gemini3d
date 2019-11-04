# these policies avoid deprecated behavior and avoid needless warning flags
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0075 NEW)
cmake_policy(SET CMP0076 NEW)

#=== compiler setup
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_gnu.cmake)
endif()

#=== MPI
find_package(MPI REQUIRED COMPONENTS Fortran)
