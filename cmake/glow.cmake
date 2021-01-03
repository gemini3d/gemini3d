include(FetchContent)

FetchContent_Declare(GLOW_proj
  GIT_REPOSITORY ${glow_url}
  GIT_TAG ${glow_tag}
  GIT_SHALLOW true
)

FetchContent_MakeAvailable(GLOW_proj)
