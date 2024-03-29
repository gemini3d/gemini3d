# CMake Install

In general on Linux, MacOS, and Windows the easiest CMake install is:

```sh
pip install cmake
```

On **macOS** the latest cmake can be obtained from the package manager:

```sh
brew install cmake
```

Most **Linux** distributions (even recent ones) come with an outdated cmake package.
There are several easy ways to get an updated CMake:

```sh
snap install cmake
```

OR, install the distribution cmake package (which will be outdated) and use it with our update script to get a newer version, e.g.:

```sh
apt install cmake
git clone https://github.com/gemini3d/external
cmake -P external/scripts/install_cmake.cmake
```

OR, install cmake from source by first downloading and extracting the latest source from [cmake.org](https://cmake.org):

```sh
cd <cmake source code directory>
./bootstrap --prefix=<install directory>
make -j
make install
```
