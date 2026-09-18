# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.25.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(ffilesystem-populate NONE)



include(ExternalProject)
ExternalProject_Add(ffilesystem-populate
                     "UPDATE_DISCONNECTED" "True" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "URL" "https://github.com/ECLAIRWaveS/ffilesystem/archive/8a9ad4b213a8fbdb19816451a714f54017170109.tar.gz" "URL_HASH" "SHA256=b0d9dfa00f6fa0977209352cff03634d71999b5f8ab4274b4ce07e56d6953e5c"
                    SOURCE_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-src"
                    BINARY_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/ffilesystem-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


