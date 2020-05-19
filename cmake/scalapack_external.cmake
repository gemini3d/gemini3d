if(NOT autobuild)
  message(STATUS "NOT autobuilding Scalapack per user -Dautobuild=off")
  return()
endif()
if(NOT SCALAPACK_FOUND)
  message(STATUS "AUTOBUILD: SCALAPACK")
endif()

include(FetchContent)

FetchContent_Declare(scalapack_proj
  GIT_REPOSITORY https://github.com/scivision/scalapack.git
  GIT_TAG v2.1.0.7
  CMAKE_ARGS "-Darith=${arith}"
)

FetchContent_MakeAvailable(scalapack_proj)

set(SCALAPACK_LIBRARIES scalapack::scalapack)
set(BLACS_LIBRARIES scalapack::blacs)
set(SCALAPACK_FOUND true)
set(BLACS_FOUND true)
