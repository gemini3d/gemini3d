
list(APPEND names h5fortran hdf5 zlib glow lapack hwm14 msis2 mumps nc4fortran scalapack matgmeini pygemini)

file(READ ${CMAKE_CURRENT_LIST_DIR}/libraries.json _libj)

foreach(n ${names})
  foreach(t url git tag zip sha1 sha256)
    string(JSON m ERROR_VARIABLE e GET ${_libj} ${n} ${t})
    if(m)
      set(${n}_${t} ${m})
    endif()
  endforeach()
endforeach()
