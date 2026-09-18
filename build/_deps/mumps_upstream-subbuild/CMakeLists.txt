# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.25.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(mumps_upstream-populate NONE)



include(ExternalProject)
ExternalProject_Add(mumps_upstream-populate
                     "UPDATE_DISCONNECTED" "True" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "URL" "https://mumps-solver.org/MUMPS_5.9.1.tar.gz" "URL_HASH" "SHA256=659c9b57646b5a003ac618baa1faf9dd2044e46c732b3daaccbc7158003e1b46"
                    SOURCE_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-src"
                    BINARY_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps_upstream-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


