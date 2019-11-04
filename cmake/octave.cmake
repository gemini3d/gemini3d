if(NOT DEFINED OctaveOK)

unset(_octpath)
if(WIN32)
  if(NOT Octave_ROOT)
    set(Octave_ROOT $ENV{HOMEDRIVE}/Octave)
  endif()
  file(GLOB _octpath "${Octave_ROOT}/Octave*/mingw64/bin/")
endif()

find_program(Octave_EXECUTABLE
  NAMES octave-cli octave
  DOC "GNU Octave"
  PATHS ${Octave_ROOT}
  HINTS ${_octpath}
  PATH_SUFFIXES bin)

# https://octave.sourceforge.io/octave/function/exist.html
# check_octave_source_runs("assert(exist('validateattributes', 'file')==2)")
# we made validateattr() so Octave 3.8 can work for compare_all()
if(Octave_EXECUTABLE)
  message(STATUS "Found GNU Octave ${Octave_EXECUTABLE}")
  set(OctaveOK true CACHE BOOL "GNU Octave is present.")
else()
  set(OctaveOK false CACHE BOOL "GNU Octave is NOT present.")
endif(Octave_EXECUTABLE)

endif()