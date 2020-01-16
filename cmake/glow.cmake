include(FetchContent)

# set(FETCHCONTENT_FULLY_DISCONNECTED true BOOL "don't download")

FetchContent_Declare(GLOW_proj
  GIT_REPOSITORY https://github.com/gemini3d/glow.git
  GIT_TAG 915592c
)

FetchContent_MakeAvailable(GLOW_proj)
