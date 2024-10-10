# Developer debug build

Array bounds checking is disabled by default due to the runtime slowdown (performance impacts).
Array bounds checking can help find bugs that cause memory corruption and intermittent crashes.

"--preset" using
[CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
is an optional way to manipulate several CMake flags via CMakePresets.json file, including
[CMAKE_BUILD_TYPE](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html).

## Single configurations

By default, Gemini3D builds with full optimizations in Release mode.
Gemini3D runs much slower (like 10x or more slower) without Release optimizations.
Avoid `-Ofast` style options that break floating-point guarantees, as Gemini3D may fail to run correctly as its internal sanity checks are broken by those flags.
Changing the build mode in a CMake project (unless using Multi config) requires reconfiguring and rebuilding each time.

Debug mode:

```sh
cmake --preset debug -Bbuild

cmake --build build
```

## Multi config (fast debug/release switching)

We use the
[Ninja](https://github.com/ninja-build/ninja/releases)
[CMake generator](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html)
to speedily switch between debug and release builds under one build directory.
Ninja is installed simply by "pip install ninja", or "brew install ninja", or by downloading, extracting, and adding the *ninja executable directory* to PATH environment variable.
Ninja speeds up the build, and allows switching rapidly between "debug" and "release" modes, where "release" is highly optimized for fastest run.

```sh
cmake --preset multi
```

sets up CMake using
[Ninja Multi Config](https://cmake.org/cmake/help/latest/generator/Ninja%20Multi-Config.html)
generator, which allows switching between Debug and Release builds at build time.
CMake sets up distinct build directories "build/Debug" and "build/Release" housing the binaries, akin to ordinary CMake builds.

Specify which to build like:

```sh
cmake --build --preset debug
# or
cmake --build --preset release
```

When using the Multi Config, one also specifies so in the self tests:

```sh
ctest --preset debug
# or
ctest --preset release
```
