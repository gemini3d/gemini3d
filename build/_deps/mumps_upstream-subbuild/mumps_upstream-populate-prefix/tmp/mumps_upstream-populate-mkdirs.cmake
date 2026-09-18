# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-build"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/tmp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/src/mumps_upstream-populate-stamp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/src/mumps_upstream-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/src/mumps_upstream-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-subbuild/mumps_upstream-populate-prefix/src/mumps_upstream-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
