if(CMAKE_VERSION VERSION_LESS 3.24 AND
  (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel" OR DEFINED ENV{MKLROOT}))
    message(WARNING "Intel oneAPI / oneMKL requires CMake >= 3.24 to take effect. LAPACK will probably not work")
endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  include(${CMAKE_CURRENT_LIST_DIR}/intel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include(${CMAKE_CURRENT_LIST_DIR}/gnu.cmake)
endif()
