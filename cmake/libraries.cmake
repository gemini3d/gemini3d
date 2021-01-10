# this reads libraries.json in memory, as a "single source of truth"
if(CMAKE_VERSION VERSION_LESS 3.19)
  # FIXME: we should eventually require CMake 3.19 for this and other stability enhancements.

  message(STATUS "Due to CMake < 3.19, using fallback Gemini library versions in ${CMAKE_CURRENT_LIST_FILE}")

  set(glow_url https://github.com/gemini3d/glow.git)
  set(glow_tag v0.981.0.1)

  set(lapack_url https://github.com/scivision/lapack.git)
  set(lapack_tag v3.9.0.2)

  set(h5fortran_url https://github.com/geospace-code/h5fortran.git)
  set(h5fortran_tag v3.4.5)

  set(mumps_url https://github.com/scivision/mumps.git)
  set(mumps_tag v5.3.5.2)

  set(nc4fortran_url https://github.com/geospace-code/nc4fortran.git)
  set(nc4fortran_tag v1.1.2)

  set(scalapack_url https://github.com/scivision/scalapack.git)
  set(scalapack_tag v2.1.0.11)

  set(pygemini_url https://github.com/gemini3d/pygemini.git)
  set(pygemini_tag main)

  return()
endif()

# preferred method CMake >= 3.19
file(READ ${CMAKE_CURRENT_LIST_DIR}/libraries.json _libj)

string(JSON glow_url GET ${_libj} "glow" "url")
string(JSON glow_tag GET ${_libj} "glow" "tag")

string(JSON lapack_url GET ${_libj} "lapack" "url")
string(JSON lapack_tag GET ${_libj} "lapack" "tag")

string(JSON h5fortran_url GET ${_libj} "h5fortran" "url")
string(JSON h5fortran_tag GET ${_libj} "h5fortran" "tag")

string(JSON mumps_url GET ${_libj} "mumps" "url")
string(JSON mumps_tag GET ${_libj} "mumps" "tag")

string(JSON nc4fortran_url GET ${_libj} "nc4fortran" "url")
string(JSON nc4fortran_tag GET ${_libj} "nc4fortran" "tag")

string(JSON scalapack_url GET ${_libj} "scalapack" "url")
string(JSON scalapack_tag GET ${_libj} "scalapack" "tag")

string(JSON pygemini_url GET ${_libj} "pygemini" "url")
string(JSON pygemini_tag GET ${_libj} "pygemini" "tag")
