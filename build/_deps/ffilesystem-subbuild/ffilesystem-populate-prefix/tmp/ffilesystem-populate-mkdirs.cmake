# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-build"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/tmp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/src/ffilesystem-populate-stamp"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/src"
  "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/src/ffilesystem-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/src/ffilesystem-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-subbuild/ffilesystem-populate-prefix/src/ffilesystem-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
