include(FetchContent)

FetchContent_Declare(scalapack_proj
  GIT_REPOSITORY https://github.com/scivision/scalapack.git
  GIT_TAG v2.1.0.3
  CMAKE_ARGS "-Darith=${arith}"
)

FetchContent_MakeAvailable(scalapack_proj)

if(NOT SCALAPACK_FOUND)
  message(STATUS "AUTOBUILD: SCALAPACK")
endif()

set(SCALAPACK_LIBRARIES scalapack::scalapack)
set(BLACS_LIBRARIES scalapack::blacs)
set(SCALAPACK_FOUND true)
set(BLACS_FOUND true)
