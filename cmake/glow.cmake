# Leave GLOW as FetchContent as we use target properties to pass up DATADIR to src/ionization
include(FetchContent)

FetchContent_Declare(GLOW_proj
  GIT_REPOSITORY ${glow_url}
  GIT_TAG ${glow_tag}
  GIT_SHALLOW true
)

FetchContent_MakeAvailable(GLOW_proj)
