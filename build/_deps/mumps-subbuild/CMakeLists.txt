# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.25.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(mumps-populate NONE)



include(ExternalProject)
ExternalProject_Add(mumps-populate
                     "UPDATE_DISCONNECTED" "True" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "URL" "https://github.com/scivision/mumps-superbuild/archive/31999e0f398d7f486c23c265ba5c40f8e1288a3c.tar.gz" "URL_HASH" "SHA256=a4cc5c87d618543af9e41ef1bea877b08bb4b3160e36bf2c602110fe4519aea2"
                    SOURCE_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-src"
                    BINARY_DIR          "/home/runner/work/gemini3d/gemini3d/build/_deps/mumps-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


