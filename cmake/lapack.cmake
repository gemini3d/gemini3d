# Finds Lapack, tests, and if not found or broken, autobuild Lapack

find_package(LAPACK)

set(lapack_external false)
if(NOT LAPACK_FOUND)
  include(${CMAKE_CURRENT_LIST_DIR}/lapack_external.cmake)
  set(lapack_external true)
endif()


if(lapack_external)
# can't run prebuild test with external libraries not yet built.
  return()
endif()
# -- verify Lapack links

set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})

if("d" IN_LIST arith)
  set(_code "double precision, external :: disnan; print *, disnan(0.); end")
elseif("s" IN_LIST arith)
  set(_code "real, external :: sisnan; print *, sisnan(0.); end")
else()
  message(WARNING "SKIP: Lapack test unknown arith: ${arith}")
  return()
endif()

check_fortran_source_compiles("${_code}" LAPACK_Compiles_OK SRC_EXT f90)
if(NOT LAPACK_Compiles_OK)
  message(WARNING "Lapack ${LAPACK_LIBRARIES} not building with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
  include(${CMAKE_CURRENT_LIST_DIR}/lapack_external.cmake)
  set(lapack_external true)
endif()