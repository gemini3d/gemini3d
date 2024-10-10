# Build Gemini3D with GCC on Linux

This method also works for
[Windows WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

```sh
apt install cmake
# or
dnf install cmake
```

If CMake is too old, install a [new CMake](./Readme_cmake_install.md).

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```

## Troubleshooting

The compiler relies on libc and libstdc++.
If build errors about "filesystem" at C++ link time, the system configuration may be messed up.
If problems, ensure environment variable LD_LIBRARY_PATH has first the libc and libstdc++ for the GCC version you wish to use.
