# Developer debug build

Those adding or modifying features inside Gemini3D code may desire to use the debugging features such as array bound checking that are disabled by default.
These debugging features make Gemini3D runs take significantly longer, but may help Gemini3D developers uncover problems with modified or added code.

If one wishes to frequently switch between Release and Debug builds, we suggest the "Multi config" section below.
Otherwise, switching from Debug to Release or vice versa requires reconfiguring the CMake project each time, which can be tedious if done frequently--and you may forget which mode you're in.

Note: the "--preset" option is using
[CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
to manipulate several CMake flags via CMakePresets.json file, including
[CMAKE_BUILD_TYPE](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html).

## Single config (must reconfigure manually to switch)

By default, Gemini3D builds with full optimizations in Release mode.
NOTE: In general one should not use `-Ofast` style options that break floating-point guarantees, as Gemini3D may fail to run correctly as its internal sanity checks are broken by those flags.
Changing the build mode in a CMake project (unless using Multi config) requires reconfiguring and rebuilding each time.

Debug mode:

```sh
cmake -B build --preset debug

cmake --build build
```

Release mode:

```sh
cmake -B build --preset release

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
cmake -B --preset multi
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
