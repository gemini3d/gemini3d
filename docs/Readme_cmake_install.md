# CMake Install

In general on Linux, MacOS, and Windows the easiest CMake install is:

```sh
pip install cmake
```

With **Homebrew** the latest cmake is available:

```sh
brew install cmake
```

**Linux** distributions generally have CMake package, available like:

```sh
snap install cmake
# or
apt install cmake
# or
dnf install cmake
```

---

It's rarely necessary to **build cmake from source**.
by first downloading and extracting the latest source from [cmake.org](https://cmake.org):

```sh
mkdir /tmp/cmake
cd /tmp/cmake

curl -LO <url of cmake source>
tar xf <cmake archive>

./bootstrap --prefix=<install directory>
make -j
make install
```
