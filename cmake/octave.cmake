function(check_octave_source_runs code)

if(NOT Octave_EXECUTABLE)
  set(ok false)
else()
  execute_process(COMMAND ${Octave_EXECUTABLE} --eval ${code}
    ERROR_QUIET OUTPUT_QUIET
    RESULT_VARIABLE ret
    TIMEOUT 10)
  if(ret EQUAL 0)
    set(OctaveOK true CACHE BOOL "GNU Octave is sufficiently new to run self-tests")
  else()
    set(OctaveOK false CACHE BOOL "GNU Octave is NOT sufficiently new to run self-tests")
  endif()
endif()
endfunction(check_octave_source_runs)


function(find_octave)
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
endif()
endfunction()

find_octave()