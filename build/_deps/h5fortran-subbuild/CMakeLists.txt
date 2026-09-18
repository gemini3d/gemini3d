# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.25.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(h5fortran-populate NONE)



include(ExternalProject)
ExternalProject_Add(h5fortran-populate
                     "UPDATE_DISCONNECTED" "True" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "URL" "https://github.com/ECLAIRWaveS/h5fortran/archive/04c7dd39405b43eac922319b6f4b74ea04b97043.tar.gz" "URL_HASH" "SHA256=9b1f546c65690ba782e7973f10b7d14b098d5f9bab83c90655a7b7e151f52281"
                    SOURCE_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-src"
                    BINARY_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/h5fortran-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


