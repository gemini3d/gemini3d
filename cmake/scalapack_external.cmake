set(scalapack_external true CACHE BOOL "autobuild ScaLapack")

include(FetchContent)

FetchContent_Declare(scalapack_proj
  GIT_REPOSITORY https://github.com/scivision/scalapack.git
  GIT_TAG v2.1.0.9
  CMAKE_ARGS "-Darith=${arith}"
)

FetchContent_MakeAvailable(scalapack_proj)

set(SCALAPACK_FOUND true)
set(BLACS_FOUND true)
