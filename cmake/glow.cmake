include(ExternalProject)

ExternalProject_Add(GLOW_proj
  GIT_REPOSITORY https://github.com/gemini3d/glow.git
  GIT_TAG master
  INSTALL_COMMAND ""  # disables the install step for the external project
)

ExternalProject_Get_Property(GLOW_proj BINARY_DIR)
set(GLOW_BINARY_DIR ${BINARY_DIR} CACHE PATH "path to GLOW")
