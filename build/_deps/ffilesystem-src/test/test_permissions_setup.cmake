if(DEFINED perm_noread)

  file(REMOVE ${perm_noread})
      # need this logic because file is non-readable and
      # kwSys:SystemTools:FileExists is false for non-readable files
  file(TOUCH ${perm_noread})
  file(CHMOD ${perm_noread} FILE_PERMISSIONS OWNER_EXECUTE)

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
    if(IS_READABLE ${perm_noread})
      message(WARNING "${perm_noread} should not be readable")
      cmake_language(EXIT 77)
    endif()
  endif()
endif()

if(DEFINED perm_nowrite)

  file(REMOVE ${perm_nowrite})
  file(TOUCH ${perm_nowrite})
  file(CHMOD ${perm_nowrite} FILE_PERMISSIONS OWNER_READ)

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
    if(IS_WRITABLE ${perm_nowrite})
      message(WARNING "${perm_nowrite} should not be writable")
      cmake_language(EXIT 77)
    endif()
  endif()
endif()


if(DEFINED perm_noexe)

  file(REMOVE ${perm_noexe})
      # need this logic because file is non-readable and
      # kwSys:SystemTools:FileExists is false for non-readable files
  file(TOUCH ${perm_noexe})
  file(CHMOD ${perm_noexe} FILE_PERMISSIONS OWNER_WRITE)

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
    if(IS_EXECUTABLE ${perm_noexe})
      message(WARNING "${perm_noexe} should not be executable")
      cmake_language(EXIT 77)
    endif()
  endif()
endif()
