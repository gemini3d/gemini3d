# MATLAB Ffilesystem FFI interface

This technique shows a general approach to using string input/output with loadlibrary() FFI for MATLAB with Ffilesystem.
It is vitally important that the compiler used to build Ffilesystem is
[MATLAB compatible for the specific MATLAB Release](https://www.mathworks.com/support/requirements/supported-compilers-linux.html).
Otherwise segmentation faults may occur when calling Ffilesystem functions from MATLAB.
This approach is shaky enough that it shouldn't be used as a default solution.
That's why in
[stdlib for Matlab](https://github.com/geospace-code/matlab-stdlib)
we use non-compiled backends like Java or .NET where native Matlab code isn't capable of doing advanced system functions.

First, build Ffilesystem shared library with CMake

```sh
cmake --workflow matlab
```

Then, use code like in
[matlab_filesystem.m](./matlab_filesystem.m)
to load the
[Ffilesystem class](./Ffilesystem.m)
and call functions.
