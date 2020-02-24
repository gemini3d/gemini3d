include(FetchContent)

FetchContent_Declare(nc4fortran_proj
  GIT_REPOSITORY https://github.com/scivision/nc4fortran.git
  GIT_TAG v0.2.0
)

FetchContent_MakeAvailable(nc4fortran_proj)
