include(FetchContent)

FetchContent_Declare(GLOW_proj
  GIT_REPOSITORY ${gemini_glow_url}
  GIT_TAG ${gemini_glow_tag}
  GIT_SHALLOW true
)

FetchContent_MakeAvailable(GLOW_proj)
