# Leave GLOW as FetchContent as we use target properties to pass up DATADIR to src/ionization
include(FetchContent)

FetchContent_Declare(GLOW
  GIT_REPOSITORY ${glow_git}
  GIT_TAG ${glow_tag})
FetchContent_MakeAvailable(GLOW)
