if(NOT package)
  return()
endif()

# make empty ref JSON file
# we will read this file and overwrite with new JSON for each sim
cmake_path(APPEND arc_dir ${PROJECT_BINARY_DIR} "upload")
cmake_path(APPEND ref_json_file ${arc_dir} "ref_data.json")

if(NOT EXISTS ${ref_json_file})
  # make a blank JSON file
  file(MAKE_DIRECTORY ${arc_dir})
  file(WRITE ${ref_json_file} "{}")
endif()

# record system metadata
file(READ ${ref_json_file} ref_json)

# check if tag exists, create if not
string(JSON m ERROR_VARIABLE e GET ${ref_json} system)
if(NOT m)
  string(JSON ref_json SET ${ref_json} system "{}")
endif()

string(JSON ref_json SET ${ref_json} system cmake_version \"${CMAKE_VERSION}\")
string(JSON ref_json SET ${ref_json} system cmake_build_type \"${CMAKE_BUILD_TYPE}\")
string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} system operating_system \"${CMAKE_HOST_SYSTEM_NAME}\")
string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} system cpu \"${host_cpu}\")
string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} system memory_ram_MB ${host_ramMB})
string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} system fortran_compiler \"${CMAKE_Fortran_COMPILER_ID}:${CMAKE_Fortran_COMPILER_VERSION}\")
string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} system c_compiler \"${CMAKE_C_COMPILER_ID}:${CMAKE_C_COMPILER_VERSION}\")

# check if tag exists, create if not
string(JSON m ERROR_VARIABLE e GET ${ref_json} gemini3d)
if(NOT m)
  string(JSON ref_json SET ${ref_json} gemini3d "{}")
endif()

string(JSON ref_json SET ${ref_json} gemini3d version \"${git_rev}\")
string(JSON ref_json SET ${ref_json} gemini3d git_branch \"${git_branch}\")
string(JSON ref_json SET ${ref_json} gemini3d git_porcelain ${git_porcelain})

# check if tag exists, create if not
string(JSON m ERROR_VARIABLE e GET ${ref_json} library)
if(NOT m)
  string(JSON ref_json SET ${ref_json} library "{}")
endif()

foreach(n LAPACK SCALAPACK MUMPS HDF5 NetCDF MPI)
  if(${n}_LIBRARIES)
    string(REPLACE ";" "," l "${${n}_LIBRARIES}")
    string(TOLOWER ${n} nl)
    string(JSON ref_json ERROR_VARIABLE e SET ${ref_json} library ${nl} \"${l}\")
  endif()

endforeach()

file(WRITE ${ref_json_file} ${ref_json})
