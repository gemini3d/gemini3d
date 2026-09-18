# Ffilesystem for Android

Using CMake and the Android
[Native Development Kit (NDK)](https://github.com/android/ndk),
Ffilesystem can be built for Android from the command line.
Android NDK supports multiple operating systems and is
[installed via Android Studio](https://developer.android.com/studio/projects/install-ndk#default-version)
or Homebrew:

```sh
brew install android-ndk
```

## CTest with Android NDK

Because Android binaries cannot run on the host, the standard
`ctest` (which executes the test binary at build time to
enumerate tests) does not work when cross-compiling. Three modes are
supported:

### Mode 1 — Build only (default)

Tests are not built. Use this for quick compile/link validation:

```sh
cmake --workflow android
```

### Mode 2 — Build and register tests, no runner

Tests are compiled and registered. `ctest` will show them as registered
but cannot run them without a device or emulator.

```sh
cmake --workflow android-tests
```

or explicitly:

```sh
cmake -B build-android --toolchain android.cmake \
  -Dffilesystem_fortran=off \
  -Dffilesystem_BUILD_TESTING=on

cmake --build build-android

ctest --test-dir build-android --show-only
```

### Mode 3 — Build, discover, and run tests with a runner

If you have an Android emulator or device configured as a runner (e.g.
via `CMAKE_CROSSCOMPILING_EMULATOR` pointing to an `adb shell` wrapper) by preset.
Before using this script / tests, the Android SDK emulator must be running and connected via adb
e.g. "adb devices" should show a device.

Assuming Android SDK is installed on Linux or macOS, in a separate terminal:

1. "source test/android.sh" from this directory
2. list available emulators with "emulator -list-avds"
3. start the emulator (e.g. "emulator -avd Medium_Phone_API_35")
4. "adb devices" should show a device
5. "cmake --workflow android-tests" to build tests for Android

```sh
cmake --workflow android-tests
```
