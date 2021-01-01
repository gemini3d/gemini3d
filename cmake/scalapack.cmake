# Finds Scalapack, tests, and if not found or broken, autobuild scalapack

if(autobuild)
  find_package(SCALAPACK)
else()
  find_package(SCALAPACK REQUIRED)
endif()

if(NOT SCALAPACK_FOUND)
  set(scalapack_external true CACHE BOOL "autobuild ScaLapack")

  include(FetchContent)

  FetchContent_Declare(scalapack_proj
    GIT_REPOSITORY ${gemini_scalapack_url}
    GIT_TAG ${gemini_scalapack_tag}
    GIT_SHALLOW true
    CMAKE_ARGS "-Darith=${arith}"
  )

  FetchContent_MakeAvailable(scalapack_proj)
endif()

target_link_libraries(SCALAPACK::SCALAPACK INTERFACE LAPACK::LAPACK)
