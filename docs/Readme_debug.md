# Developer debug build

Those adding or modifying features inside Gemini3D code may desire to use the debugging features such as array bound checking that are disabled by default.
These debugging features make Gemini3D runs take significantly longer, but may help Gemini3D developers uncover problems with modified or added code.

If one wishes to frequently switch between Release and Debug builds, we suggest the "Multi config" section below.
Otherwise, switching from Debug to Release or vice versa requires reconfiguring the CMake project each time, which can be tedious if done frequently--and you may forget which mode you're in.

## Single config (must reconfigure manually to switch)

By default, Gemini3D builds with full optimizations in Release mode.
NOTE: In general one should not use `-Ofast` style options that break floating-point guarantees, as Gemini3D may fail to run correctly as its internal sanity checks are broken by those flags.
Changing the build mode in a CMake project (unless using Multi config) requires reconfiguring and rebuilding each time.

Debug mode:

```sh
cmake --preset debug

cmake --build build
```

Release mode:

```sh
cmake --preset default

cmake --build build
```

## Multi config (fast debug/release switching)

In general, we recommend developers use the [Ninja backend](https://github.com/ninja-build/ninja/releases), which may be installed simply by "pip install ninja", or "brew install ninja", or by downloading, extracting, and adding the *ninja executable directory* to PATH environment variable.
Ninja offers a significant build speedup, and allows switching rapidly between "debug" and "release" modes, where "release" is highly optimized for fastest run.

```sh
cmake --preset multi
```

sets up CMake in "Ninja Multi Config" mode, which allows switching between Debug and Release at build time.
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
