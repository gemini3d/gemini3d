set(_names glow hwm14 msis2
lapack mumps scalapack
matgmeini pygemini
nc4fortran h5fortran hdf5)

file(READ ${CMAKE_CURRENT_LIST_DIR}/libraries.json _libj)

foreach(n ${_names})
  foreach(t url git tag zip sha256)
    string(JSON m ERROR_VARIABLE e GET ${_libj} ${n} ${t})
    if(m)
      set(${n}_${t} ${m})
    endif()
  endforeach()
endforeach()

# Zlib special case to allow fallback
if(zlib_legacy)
  string(JSON zlib_url GET ${_libj} zlib1 url)
  string(JSON zlib_sha256 GET ${_libj} zlib1 sha256)
else()
  string(JSON zlib_url GET ${_libj} zlib2 url)
  string(JSON zlib_sha256 GET ${_libj} zlib2 sha256)
endif()
