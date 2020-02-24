include(FetchContent)

FetchContent_Declare(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG v2.6.2
)

FetchContent_MakeAvailable(h5fortran_proj)
