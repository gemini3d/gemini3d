# Distributed under the OSI-approved BSD 3-Clause License.
# Copyright 2013, Julien Schueller
# Copyright 2018, Michael Hirsch, Ph.D.

#[=======================================================================[.rst:
FindOctave
----------

Finds GNU Octave and provides Octave tools, libraries and compilers to CMake.

This packages primary purposes are

* find the Octave exectuable to be able to run unit tests on ``.m`` code
* to run specific commands in Octave
* to retrieve various information from Octave (versions)

Module Input Variables
^^^^^^^^^^^^^^^^^^^^^^

Users or projects may set the following variables to configure the module
behaviour:

:variable:`Octave_ROOT_DIR`
  the root of the Octave installation.

Variables defined by the module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Result variables
""""""""""""""""

``Octave_FOUND``
  ``TRUE`` if the Octave installation is found, ``FALSE`` otherwise.
``Octave_EXECUTABLE``
  Octave interpreter
``Octave_INCLUDE_DIRS``
  include path for mex.h
``Octave_LIBRARIES``
  octinterp, octave, cruft
``Octave_OCTINTERP_LIBRARY``
  path to the library octinterp
``Octave_OCTAVE_LIBRARY``
  path to the liboctave
``Octave_CRUFT_LIBRARY``
  path to the libcruft
``Octave_VERSION``
  Octave version string
``Octave_MAJOR_VERSION``
  major version
``Octave_MINOR_VERSION``
  minor version
``Octave_PATCH_VERSION``
  patch version
``Octave_OCT_FILE_DIR``
  object files that will be dynamically loaded
``Octave_OCT_LIB_DIR``
  oct libraries
``Octave_ROOT_DIR``
  octave prefix
#]=======================================================================]


find_program(Octave_CONFIG_EXECUTABLE
             NAMES octave-config
             HINTS ${Octave_ROOT_DIR})

if(Octave_CONFIG_EXECUTABLE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTAVE_HOME
                  OUTPUT_VARIABLE Octave_ROOT_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p BINDIR
                  OUTPUT_VARIABLE Octave_BIN_PATHS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTINCLUDEDIR
                  OUTPUT_VARIABLE Octave_INCLUDE_PATHS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTLIBDIR
                  OUTPUT_VARIABLE Octave_LIBRARIES_PATHS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTFILEDIR
                  OUTPUT_VARIABLE Octave_OCT_FILE_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTLIBDIR
                  OUTPUT_VARIABLE Octave_OCT_LIB_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  find_program(Octave_EXECUTABLE
               NAMES octave
               HINTS ${Octave_BIN_PATHS}
              )

  find_library(Octave_OCTINTERP_LIBRARY
             NAMES octinterp liboctinterp
             HINTS ${Octave_LIBRARIES_PATHS}
            )
  find_library(Octave_OCTAVE_LIBRARY
               NAMES octave liboctave
               HINTS ${Octave_LIBRARIES_PATHS}
              )
  find_library(Octave_CRUFT_LIBRARY
               NAMES cruft libcruft
               HINTS ${Octave_LIBRARIES_PATHS}
              )

  set(Octave_LIBRARIES ${Octave_OCTINTERP_LIBRARY})
  list(APPEND Octave_LIBRARIES ${Octave_OCTAVE_LIBRARY})
  if (Octave_CRUFT_LIBRARY)
    list(APPEND Octave_LIBRARIES ${Octave_CRUFT_LIBRARY})
  endif()

  find_path(Octave_INCLUDE_DIR
            NAMES mex.h
            HINTS ${Octave_INCLUDE_PATHS}
            )

  set(Octave_INCLUDE_DIRS ${Octave_INCLUDE_DIR})
else()
  find_program(Octave_EXECUTABLE
               NAMES octave
               HINTS ${Octave_ROOT_DIR})
endif()



execute_process(COMMAND ${Octave_EXECUTABLE} -v
                OUTPUT_VARIABLE Octave_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)

if(Octave_VERSION)
  string(REGEX REPLACE "GNU Octave, version ([0-9])\\.[0-9]+\\.[0-9]+.*" "\\1" Octave_MAJOR_VERSION ${Octave_VERSION})
  string(REGEX REPLACE "GNU Octave, version [0-9]\\.([0-9]+)\\.[0-9]+.*" "\\1" Octave_MINOR_VERSION ${Octave_VERSION})
  string(REGEX REPLACE "GNU Octave, version [0-9]\\.[0-9]+\\.([0-9]+).*" "\\1" Octave_PATCH_VERSION ${Octave_VERSION})

  set(Octave_VERSION ${Octave_MAJOR_VERSION}.${Octave_MINOR_VERSION}.${Octave_PATCH_VERSION})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Octave
  REQUIRED_VARS Octave_EXECUTABLE
  VERSION_VAR Octave_VERSION
  HANDLE_COMPONENTS)

mark_as_advanced(
  Octave_OCT_FILE_DIR
  Octave_OCT_LIB_DIR
  Octave_OCTINTERP_LIBRARY
  Octave_OCTAVE_LIBRARY
  Octave_CRUFT_LIBRARY
  Octave_LIBRARIES
  Octave_INCLUDE_DIR
  Octave_INCLUDE_DIRS
  Octave_ROOT_DIR
  Octave_VERSION
  Octave_MAJOR_VERSION
  Octave_MINOR_VERSION
  Octave_PATCH_VERSION
)
