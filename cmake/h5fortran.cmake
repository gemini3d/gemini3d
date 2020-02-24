include(FetchContent)

FetchContent_Declare(h5fortran_proj
  GIT_REPOSITORY https://github.com/scivision/h5fortran.git
  GIT_TAG v2.6.3
)

FetchContent_MakeAvailable(h5fortran_proj)

if(NOT HDF5OK)
  message(FATAL_ERROR "HDF5 was requested but it's not working on your system")
endif()