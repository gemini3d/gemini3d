# Cosmopolitan universal binary

[Cosmopolitan](https://github.com/jart/cosmopolitan)
provides a universial binary that can run on many computing platforms (CPU archs and operating systems).

To build with Cosmo, assuming the
[Cosmo binaries](https://github.com/jart/cosmopolitan?tab=readme-ov-file#getting-started)
are under $HOME/cosmo/, make a source file $HOME/cosmo.sh like:

```sh
export COSMO_ROOT=$HOME/cosmo

b=$COSMO_ROOT/bin

export CFLAGS=-I$COSMO_ROOT/include
export CXXFLAGS=-I$COSMO_ROOT/include
export LDFLAGS=-L$COSMO_ROOT/lib
export PKG_CONFIG="$b/pkg-config --with-path=$COSMO_ROOT/lib/pkgconfig"
export INSTALL=$b/cosmoinstall

export AR=$b/cosmoar
export RANLIB=$b/cosmoranlib

export ASM=$b/cosmocc
export CC=$b/cosmocc
export CXX=$b/cosmoc++
```

if you have trouble, try the CPU-specific names like "aarch64-unknown-cosmo-cc" and so on.

Then to build a universal binary with Cosmopolitan, use the following commands:

```sh
cmake --workflow cosmo
```

Run the universal binary

```sh
build-cosmo/app/fs_cli
```

The self-tests all error like "Unknown system error -8" from CTest, but can be run correctly individually from the command line.

```sh
ctest --test-dir build-cosmo
```
