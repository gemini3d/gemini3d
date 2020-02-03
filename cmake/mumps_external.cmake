include(FetchContent)

FetchContent_Declare(MUMPS_proj
  GIT_REPOSITORY https://github.com/scivision/mumps.git
  GIT_TAG v5.2.1.8
  CMAKE_ARGS "-Darith=${arith}" "-Dparallel=true" "-Dmetis=${metis}" "-Dscotch=${scotch}" "-Dopenmp=false"
)

FetchContent_MakeAvailable(MUMPS_proj)

if(NOT MUMPS_FOUND)
  message(STATUS "AUTOBUILD: MUMPS")
endif()

set(MUMPS_LIBRARIES mumps::mumps)
set(MUMPS_FOUND true)
