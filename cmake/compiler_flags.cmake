# === compiler setup
# feel free to add more compiler_*.cmake
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(WIN32)
    include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel_win.cmake)
  else(WIN32)
    include(${CMAKE_CURRENT_LIST_DIR}/compiler_intel.cmake)
  endif(WIN32)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_gnu.cmake)
endif()
