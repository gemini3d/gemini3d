# Ffilesystem Benchmarks

These were run ~ commit 8424937f using `CMAKE_BUILD_TYPE=Release`.

* which(): wallclock time is proportional to the number of directories in PATH and which position in PATH the executable happens to be at.
* homedir(): wallclock time is proportional to the number of environment variables.
* expanduser(): wallclock time varies considerably depending if "~" or "~/" are the leading characters or not. Also calls homedir().
* canonical(): wallclock time depends on if the path is already absolute and is proportional to the number of symlinks in the path, length of the path, and also calls expanduser().


## Linux

Intel Xeon workstation

### Clang AMD Clang 16.0.3 (CLANG: AOCC_4.2.0-Build#89 2023_12_13)

```
Cpp: 1000 x canonical(~/..) = /home: 2.651 us
Cpp: 1000 x expanduser(~/..) = /home: 1.201 us
Cpp: 1000 x homedir() = /home/xxx: 0.425 us
Cpp: 1000 x normal(~/..) = .: 0.36 us
Cpp: 1000 x resolve(~/..) = /home: 2.666 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 41.484 us
```

### NVIDIA nvc 24.3.0

```
Cpp: 1000 x canonical(~/..) = /home: 2.662 us
Cpp: 1000 x expanduser(~/..) = /home: 1.268 us
Cpp: 1000 x homedir() = /home/xxx: 0.422 us
Cpp: 1000 x normal(~/..) = .: 0.365 us
Cpp: 1000 x resolve(~/..) = /home: 2.827 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 13.265 us
```

### GNU GCC 13.2.1

```
Cpp: 1000 x canonical(~/..) = /home: 2.708 us
Cpp: 1000 x expanduser(~/..) = /home: 1.171 us
Cpp: 1000 x homedir() = /home/xxx: 0.415 us
Cpp: 1000 x normal(~/..) = .: 0.293 us
Cpp: 1000 x resolve(~/..) = /home: 2.746 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 15.868 us
```

### Intel LLVM 20240000 Intel(R) oneAPI DPC++/C++ Compiler 2024.0.0 (2024.0.0.20231017)

```
Cpp: 1000 x canonical(~/..) = /home: 2.699 us
Cpp: 1000 x expanduser(~/..) = /home: 1.211 us
Cpp: 1000 x homedir() = /home/xxx: 0.427 us
Cpp: 1000 x normal(~/..) = .: 0.363 us
Cpp: 1000 x resolve(~/..) = /home: 2.738 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 34.633 us
```

## Windows

Core i5 laptop

### GNU GCC 13.2.0

```
Cpp: 1000 x canonical(~/..) = C:/Users: 94.2 us
Cpp: 1000 x expanduser(~/..) = C:/Users: 4.6 us
Cpp: 1000 x homedir() = C:/Users/xxx: 2.5 us
Cpp: 1000 x normal(~/..) = .: 0.4 us
Cpp: 1000 x resolve(~/..) = C:/Users: 93.2 us
Cpp: 1000 x which(cmd.exe) = C:/WINDOWS/system32/cmd.exe: 301.2 us
```

### MSVC 193933522

```
Cpp: 1000 x canonical(~/..) = C:/Users: 37.9 us
Cpp: 1000 x expanduser(~/..) = C:/Users: 4.4 us
Cpp: 1000 x homedir() = C:/Users/xxx: 3.3 us
Cpp: 1000 x normal(~/..) = .: 0.2 us
Cpp: 1000 x resolve(~/..) = C:/Users: 38 us
Cpp: 1000 x which(cmd.exe) = C:/WINDOWS/system32/cmd.exe: 160.7 us
```

### Intel LLVM 20240002 Intel(R) oneAPI DPC++/C++ Compiler 2024.0.2 (2024.0.2.20231213)

```
Cpp: 1000 x canonical(~/..) = C:/Users: 40.3 us
Cpp: 1000 x expanduser(~/..) = C:/Users: 5.7 us
Cpp: 1000 x homedir() = C:/Users/xxx: 4.7 us
Cpp: 1000 x normal(~/..) = .: 0.2 us
Cpp: 1000 x resolve(~/..) = C:/Users: 40.3 us
Cpp: 1000 x which(cmd.exe) = C:/WINDOWS/system32/cmd.exe: 272.7 us
```

### WSL (Ubuntu 22.04): GNU GCC 11.4.0

```
Cpp: 1000 x canonical(~/..) = /home: 2.116 us
Cpp: 1000 x expanduser(~/..) = /home: 0.897 us
Cpp: 1000 x homedir() = /home/xxx: 0.304 us
Cpp: 1000 x normal(~/..) = .: 0.212 us
Cpp: 1000 x resolve(~/..) = /home: 2.129 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 3.213 us
```

### Cygwin: GNU GCC 11.4.0

```
Cpp: 1000 x canonical(~/..) = /home: 105.9 us
Cpp: 1000 x expanduser(~/..) = /home: 4.4 us
Cpp: 1000 x homedir() = /home/xxx: 1.6 us
Cpp: 1000 x normal(~/..) = .: 0.8 us
Cpp: 1000 x resolve(~/..) = /home: 100.3 us
Cpp: 1000 x which(sh) = /usr/bin/sh: 92.2 us
```
