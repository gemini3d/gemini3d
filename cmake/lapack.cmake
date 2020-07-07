# Finds Lapack, tests, and if not found or broken, autobuild Lapack
if(lapack_external)
  include(${CMAKE_CURRENT_LIST_DIR}/lapack_external.cmake)
  return()
endif()

find_package(LAPACK)

if(NOT LAPACK_FOUND)
  include(${CMAKE_CURRENT_LIST_DIR}/lapack_external.cmake)
  set(lapack_external true CACHE BOOL "autobuild Lapack")
  return()
else()
  set(lapack_external false CACHE BOOL "autobuild Lapack")
endif()

# -- verify Lapack links

set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)

if("d" IN_LIST arith)
  set(_code "double precision, external :: disnan; print *, disnan(0.); end")
elseif("s" IN_LIST arith)
  set(_code "real, external :: sisnan; print *, sisnan(0.); end")
else()
  message(STATUS "SKIP: Lapack test arith: ${arith}")
  return()
endif()

include(CheckFortranSourceCompiles)
check_fortran_source_compiles("${_code}" LAPACK_Compiles_OK SRC_EXT f90)
if(NOT LAPACK_Compiles_OK)
  message(STATUS "Lapack ${LAPACK_LIBRARIES} not building with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  include(${CMAKE_CURRENT_LIST_DIR}/lapack_external.cmake)
  set(lapack_external true CACHE BOOL "autobuild Lapack")
endif()
