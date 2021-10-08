# === compiler setup
# feel free to add more *.cmake
if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  include(${CMAKE_CURRENT_LIST_DIR}/intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)
endif()
