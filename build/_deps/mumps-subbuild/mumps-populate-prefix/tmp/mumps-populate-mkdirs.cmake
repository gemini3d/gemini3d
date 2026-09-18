# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-build"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/tmp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/src/mumps-populate-stamp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/src/mumps-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/src/mumps-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-subbuild/mumps-populate-prefix/src/mumps-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
