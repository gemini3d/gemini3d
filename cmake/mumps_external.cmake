set(mumps_external true CACHE BOOL "autobuild Mumps")

include(FetchContent)

FetchContent_Declare(MUMPS_proj
  GIT_REPOSITORY https://github.com/scivision/mumps.git
  GIT_TAG v5.3.3.5
  CMAKE_ARGS "-Darith=${arith}" "-Dparallel=true" "-Dmetis=${metis}" "-Dscotch=${scotch}" "-Dopenmp=false"
)

FetchContent_MakeAvailable(MUMPS_proj)

set(MUMPS_FOUND true)
