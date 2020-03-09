if(NOT DEFINED octave_disabled)

unset(_octpath)
if(WIN32)
  if(NOT Octave_ROOT)
    if(DEFINED ENV{Octave_ROOT})
      set(Octave_ROOT $ENV{Octave_ROOT})
    else()
      set(Octave_ROOT $ENV{HOMEDRIVE}/Octave)
    endif()
  endif()
  file(GLOB _octpath "${Octave_ROOT}/Octave*/mingw64/bin/")
endif()

find_program(Octave_EXECUTABLE
  NAMES octave-cli octave
  DOC "GNU Octave"
  PATHS ${Octave_ROOT} ENV ${OCTAVE_EXECUTABLE}
  HINTS ${_octpath}
  PATH_SUFFIXES bin)

if(Octave_EXECUTABLE)
  message(STATUS "Found GNU Octave ${Octave_EXECUTABLE}")
  set(octave_disabled false CACHE BOOL "GNU Octave tests")
else()
  set(octave_disabled true)
endif()

endif()
