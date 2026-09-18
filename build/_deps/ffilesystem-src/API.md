# Fortran Filesystem API

Fortran filesystem modules contains numerous procedures and one (optional, default enabled) Fortran type "path_t" that contains properties and methods.

C++ stdlib `<filesystem>` is used extensively within Ffilesystem to implement functions in a platform-agnostic and robust way.
Fallback to plain C++17 is available for compilers that do not support C++
[<filesystem>](https://en.cppreference.com/w/cpp/header/filesystem).
For the interchange of character strings between Fortran and C++ / C, the buffer length is determined at compile time and is available in `fs_get_max_path()` (C, C++) or `max_path()` (Fortran).

```fortran
integer :: m
m = max_path()
```

The maximum length of path segments are also limited, typically to 255 characters.
That is, the overall path might be allowed to be thousands of characters long, but each segment is limited individually as well.

```fortran
integer :: m
m = max_component("/")
```

For Windows the top-level drive "D:/" and similar requires a trailing slash, as in numerous other programs.

## path_t

`path_t` can be manually disabled with CMake by setting `cmake -DHAVE_F03TYPE=0`.

The "path_t" type uses getter and setter procedure to access the path as a string `character(:), allocatable`.

```fortran
use filesystem, only : path_t

type(path_t) :: p

p = path_t("my/path")  !< setter

print *, "path: ", p%path() !< getter
```

The retrieved path string may be indexed like:

```fortran
p%path(2,4)  !< character index 2:4

p%path(2) !< character index 2:end
```

In all the examples, we assume "p" is path_t.

## System capabilities

Character, allocatable: the Fortran compiler name and version

```
compiler()
```

---

Character, allocatable: the C/C++ compiler name and version.

```
compiler_c()
```

---

Character, allocatable: the Terminal shell that called the program

```
get_shell()
```

On Windows systems the shell is typically "cmd" or "powershell" (pwsh).
On non-Windows systems (or WSL) the shell is typically like "bash" or "zsh".

---

Character, allocatable: the Terminal emulator that called the program

```fortran
get_terminal()
```

The Terminal on Windows is "ConsoleWindowClass" for legacy COMSPEC or "PseudoConsoleWindow" for Windows Terminal.
On non-Windows systems (or WSL) the Terminal may represent itself like "xterm" or "gnome-terminal".

---

Character, allocatable: the CPU architecture

```fortran
arch = cpu_arch()
```

---

character: Ffilesystem backend (e.g. <filesystem> or C)

```fortran
fs_backend()
```

---

Ffilesystem has `ISO_Fortran_binding.h` available and using it for at least some functions

```fortran
fs_has_cfi()
```


Logical: ffilesystem was compiled with optimizations:

```fortran
fs_is_optimized()
```

Logical: path of filesystem detected as case sensitive. Path must be writable for check.

```fortran
is_case_sensitive("/my/path")
```

integer (long): C++ level of macro `__cplusplus`

```fortran
cpp_lang()
```

integer (long): C `__STDC_VERSION__` macro

```fortran
c_lang()
```

## subroutines

These subroutines are available in the "filesystem" module.

```fortran
call create_symlink("my/path", "my/symlink", ok)

logical, intent(out), optional :: ok !< true if succeeded
```

Copy source to destination.
Optionally, overwrite existing file.

```fortran
character(*) :: dest = "new/file.ext"

call p%copy_file(dest)
! or
call copy_file("original.txt", "acopy.txt")

! overwrite
call copy_file("original.txt", "acopy.txt", overwrite=.true.)

character(*), intent(in) :: source, dest
logical, intent(in), optional :: overwrite
logical, intent(out), optional :: ok !< true if successful
```

Make directory with parent directories if specified

```fortran
p = path_t("my/new/dir")
! suppose only directory "my" exists
call p%mkdir()
! now directory my/new/dir exists
! OR
call mkdir("my/new/dir")
```

Touch file (create empty file if not a file).
The directories containing the file must already exist.
Also updates the file access/modification times to current time.

```fortran
logical, optional :: ok
type(path_t) :: p

call p%touch(ok)

call touch("myfile.ext", ok)
```

Get last modified time of path.

```fortran
use, intrinsic :: iso_c_binding, only : C_LONG
integer(C_LONG) :: mtime

mtime = get_modtime("my/file.txt")
```

Set modified time of path to current time.

```fortran
logical :: ok

ok = set_modtime("my/file.txt")
```


Delete file, empty directory, or symbolic link (the target of a symbolic link is not deleted).

```fortran
call p%remove()
! or
call remove("my/file.txt")
```

Rename file or directory. Target "to" will be overwritten if it exists.

```fortran
fs_rename(from, to)
```

create symbolic link to file or directory:

```fortran
call p%create_symlink(link)
! or
call create_symlink(target, link)
```

## path_t

These methods emit a new "path_t" object.
It can be a new path_t object, or reassign to the existing path_t object.

On Windows, force file separators (if any) to POSIX "/"

```fortran
p = path_t('my\path')
p = p%as_posix()

! my/path
```

Expand home directory, swapping file separators "\" for "/" and drop redundant file separators "//".

```fortran
! Fortran does not understand tilde "~"

p = path_t("~/my/path")
p = p%expanduser()
```

Read symlink target if path is a symbolic link--empty string if not a symlink.

```fortran
target = p%read_symlink()
! or
target = read_symlink("my/symlink")
```

Realpath(): usually users will want resolve() or canonical() instead

```fortran
character(:), allocatable :: realpath(".././mypath")
```

Resolve path. This means to absolutize, canonicalize, resolving symbolic links.

* "strict" if true required the path to exist (default false).

```fortran
function resolve(path, strict)
character(*), intent(in) :: path
logical, intent(in), optional :: strict

p = path_t("~/../b")
p = p%resolve()

p%path() == "<absolute path of user home directory>/b"

! --- relative path resolved to current working directory
p = path_t("../b")
p = p%resolve()

p%path() == "<absolute path of current working directory>/b"
```

Canonicalize path. This means to normalize, resolve symbolic links, resolve relative paths.
For case-insensitive filesystems, the path's actual case is returned.
If the path doesn't exist and no absolute path is given, the path is resolved as far as possible with existing path components, and then ".", ".." are lexicographically resolved.
The non-existing path may be absolutized based on the current working directory depending on the system (macOS does this).

* "strict" if true required the path to exist (default false).

```fortran
function canonical(path, strict)
character(*), intent(in) :: path
logical, intent(in), optional :: strict

p = path_t("~/../b")
p = p%canonical()

p%path() == "<absolute path of user home directory>/b"

! --- relative path resolved to current working directory
p = path_t("../b")
p = p%canonical()

p%path() == "../b"
```

Swap file suffix

```fortran
p = path_t("my/file.h5")

p = p%with_suffix(".hdf5")

! p%path() == "my/file.hdf5"
```

Normalize path, a lexical operation removing ".." and "." and duplicate file separators "//".
The path need not exist.
Trailing file separators are gobbled.

```fortran
p = p%normal()
! or
normal("./my//path/../b/")  !< "my/b"
```

Join path with other path string using posix separators.
The paths are treated like strings.
No path resolution is used, so non-sensical paths are possible for non-sensical input.

```fortran
p = path_t("a/b")

p = p%join("c/d")

! p%path == "a/b/c/d"
```

## integer(int64)

These procedures emit a 64-bit integer value.

len_trim() of p%path()

```fortran
p%length()
```

File size (bytes):

```fortran
p%file_size()
! or
file_size("my/file.txt")
```

Maximum number of open files for this process:

```fortran
get_max_open_files()
```

Space available to user on drive containing path (bytes):

```fortran
space_available("c:/")
```

Drive capacity available to user on drive (bytes):

```fortran
space_capacity("c:/")
```

Block size of filesystem (bytes). For example, a typical block size is 4096 bytes.

```fortran
get_blksize("/")
```

Hard link count of file or directory:

```fortran
hard_link_count("my/file.txt")
```

## logical

These methods emit a logical value.


Is directory or file empty. False if path doesn't exist

```fortran
is_empty("my/dir")
```

Does directory exist:

```fortran
p%is_dir()
! or
is_dir("my/dir")
```

Does the path have a filename (lexical operation, not checking if path exists):

```fortran
has_filename("a/b/c")  !< true
has_filename("a/b/")   !< false
```


Is the drive removable (USB drive, CD/DVD drive):

```fortran
is_removable("c:/")
```

Error stop if directory does not exist

```fortran
call assert_is_dir("my/dir")
```

Is prefix of path. This is a lexical operation.

```fortran
p%is_prefix(prefix)

is_prefix("my/dir", "my")
```

Is path a subdirectory under (not just equal to) of "dir". this is a lexical operation.

```fortran
p%is_subdir(dir)

is_subdir("my/dir", "my")
```

Is filename "safe" (no path separators, no reserved names, no special characters, no white space):

```fortran
logical :: is_safe_name()

is_safe_name("my_file.txt")
```

---

Is "path" a file or directory (or a symbolic link to existing file or directory).
Like Python pathlib.Path.exists(), even if the path does not have read permission,
it still may exist.

NOTE: on Windows `exists("./not-exist/..")` is `true`, but on Unix `false`.
In C/C++ `access()` or `stat()` the same behavior is observed Windows vs Unix.

```fortran
p%exists()
! or
exists("my/file.txt")
```

---

Does path exist, or is it a symbolic link that is broken (pointing to a non-existent target)?

```fortran
lexists("my/file.txt")
```

---

is "path" a file or a symbolic link to an existing file.
Like Python pathlib.Path.is_file(), even if the file does not have read permission,
it still may exist.

```fortran
p%is_file()
! or
is_file("my/file.txt")
```

Is the path a special character device (like a terminal or /dev/null)?

```fortran
p%is_char_device()
! of
is_char_device("/dev/null")
```

Is the path a FIFO (named pipe)?

```fortran
is_fifo("my/pipe")
```

Is the path a Windows App Execution Alias?

```fortran
s1 = which("bash.exe")

is_appexec_alias(s1)
```

On Windows, is the path a reserved name (like "nul" on Windows)?

```fortran
p%is_reserved()
! or
is_reserved("nul")
```

Error stop if file does not exist

```fortran
call assert_is_file("my/dir")
```

Is path a symbolic link:

```fortran
p%is_symlink()
! or
is_symlink("my/path")
```

Is path absolute:

```fortran
p%is_absolute()
! or
is_absolute("my/path")
```

## Equivalent paths

Does path "p" *resolve* to the same path as "other".
Does not expanduser().
To be true:

* path must exist
* path must be traversable  E.g. "a/b/../c" resolves to "a/c" iff a/b also exists.
* symlink resolves to its target on the SAME drive - networked drives and Windows Dev Drives map to different drives (even if virtual drive on same physical drive), so they do not resolve to the same path and hence report "false" for same_file() or fs_equivalent() in C++.

```fortran
p%same_file(other)
! or
same_file(path1, path2)
```

## file permissions

Is file (script or binary executable) executable by the user.

* Windows: readable file with a path suffix in environment variable PATHEXT (at least, .com|.exe|.bat|.cmd) is executable.
* non-Windows: executable permission bits for the owner, group, or others file determine if it is executable.

```fortran
!! logical

p%is_exe()
! or
is_exe("my/file.exe")
```

---

Is file detected as a binary executable AND readable.

* Windows: GetBinaryType() is used to determine if a file is a binary executable.
* non-Windows: magic number of the file is checked, which is a heuristic and may not be 100% accurate.

```fortran
is_executable_binary("my/file.exe")
```

---

Make regular file executable (or not) for owner.

Windows: set_permissions(path, executable=) does NOT work (MinGW, oneAPI, MSVC).

```fortran
!! subroutine

call p%set_permissions(readable, writable, executable=.true., ok)
! or
call set_permissions("my/file.exe", executable=.true.)

logical, intent(in), optional :: readable, writable, executable
logical, intent(out), optional :: ok  !< true if successful
```

---

Is path (file or directory) readable by the user.

```fortran
!! logical

p%is_readable()
! or
is_readable("my/file.txt")
is_readable("./")
```

Is path (file or directory) writable by the user.

```fortran
!! logical

p%is_writable()
! or
is_writable("my/file.txt")
is_writable("./")
```

## character(:), allocatable

These procedures emit a string.

---

Force file separators (if any) to Posix "/"

```fortran
as_posix('my\path')
! my/path
```

Join path_t with other path string using posix separators.
The paths are treated like strings.
No path resolution is used, so non-sensical paths are possible for non-sensical input.

```fortran
join("a/b", "c/d") ! "a/b/c/d"
```

---

Cygwin-specific:

```fortran
character(:), allocatable :: to_cygpath(winpath)

character(:), allocatable :: to_winpath(cygpath)
---

transform to/from Windows path to Cygwin POSIX path

---

Find executable file on environment variable PATH if present.
Windows must include the ".exe" or ".bat" etc. suffix.
Windows prioritizes CWD.
Does not resolve path--if Windows CWD or relative path is in PATH, may output relative path.
Does not expanduser tilde.

On Windows, security / virus scanners may block cmd.exe and similar under `%SYSTEMROOT%` from being found.

```fortran
character(:), allocatable :: which("myprog")

! or

which("myprog", "/path/to/search")

! or

which("myprog", find_all=.true.)
```

`find_all` if true returns all found executables in a list separated by the system path separator.
Otherwise, the default false returns the first found executable.

---

Get file suffix: extracts path suffix, including the final "." dot

```fortran
p%suffix()
! or
suffix("my/file.txt")  !< ".txt"
```

Swap file suffix

```fortran
with_suffix("to/my.h5", ".hdf5")  !< "to/my.hdf5"
```

Get parent directory of path.
The parent of the top-most relative path is ".".
We define the parent of a path as the directory above the specified path.
Trailing slashes are gobbled.
The path is not normalized.

```fortran
p%parent()
! or
parent("my/file.txt")  !< "my"

parent("a") !< "."
```

Get file name without path:

```fortran
p%file_name()
! or
file_name("my/file.txt")  ! "file.txt"
```

Get file name without path and suffix:

```fortran
p%stem()
! or
stem("my/file.txt")  !< "file"
```

Get path separator: ":" for Unix, ";" for Windows

```fortran
character :: s = pathsep()
```

Get null path: "/dev/null" on Unix, "nul" on Windows

```fortran
character(:), allocatable :: n = devnull()
```

Get filesystem type of path.
E.g. "NTFS" or "ext4". WSL shows "v9fs" for Windows paths.

```fortran
character(:), allocatable :: p%filesystem_type()

! or

filesystem_type("my/file.txt")  !< "NTFS" or "ext4" or ...
```

Get drive root path. E.g. Unix "/"  Windows "c:/"
Requires absolute path or will return empty string.

```fortran
p%root()
! or
root("/a/b/c") !< "/" on Unix, "" on Windows

root ("c:/a/b/c") !< "c:" on Windows, "" on Unix
```

Get drive root name. On non-Windows, returns empty string.

```fortran
p%root_name()
! or
root_name("c:/a/b/c") !< "c:" on Windows, "" on Unix
```

Expand user home directory.

```fortran
expanduser("~/my/path")   !< "/home/user/my/path" on Unix, "<root>/Users/user/my/path" on Windows
```

---

Canonicalize path. This means to normalize, resolve symbolic links, and resolve relative paths when the path exists.
If the path doesn't exist and no absolute path is given, the path is resolved as far as possible with existing path components, and then ".", ".." are lexicographically resolved.

* "strict" if true required the path to exist (default false).

```fortran
function canonical(path, strict)
character(*), intent(in) :: path
logical, intent(in), optional :: strict

canonical("../b")
```

---

Resolve path.
First attempts to resolve an existing path.
If that fails, the path is resolved as far as possible with existing path components, and then ".", ".." are lexicographically resolved.

* "strict" if true required the path to exist (default false).

```fortran
function resolve(path, strict)
character(*), intent(in) :: path
logical, intent(in), optional :: strict

! --- relative path resolved to current working directory
resolve("../b")
```

## Windows-specific functions

Windows: long to short path

```fortran
shortname("C:/Program Files")  !< "C:/PROGRA~1"
```

Windows: short to long path

```fortran
longname("C:/PROGRA~1")  !< "C:/Program Files"
```

Are Windows long paths enabled (longer than 260 characters)?

```fortran
long_paths_enabled()  !< logical
```

## relative / proximate paths

Get path relative to other path, treating each path as "weakly_canonical".
This requires the `<filesystem>` backend.

```fortran
relative_to(base, other)
character(*), intent(in) :: base, other


relative_to("/a/b/c", "/a/b")  !< "c"

p = path_t("/a/b/c")
p%relative_to("/a")  !< "b/c"

p%relative_to("d")  !< ""
```

Get path proximate to other path, treating each path as "weakly_canonical".
This requires the `<filesystem>` backend.

```fortran
proximate_to(base, other)
character(*), intent(in) :: base, other

proximate_to("/a/b/c", "/a/b")  !< "c"

p = path_t("/a/b/c")
p%proximate_to("d")  !< "d"
```

## System

Get home directory, or empty string if not found.

```fortran
character(:), allocatable :: get_homedir()
```

Get profile/pw directory. `get_homedir()` is normally preferred.

```fortran
character(:), allocatable :: get_profile_dir()
```

Get hostname of computer

```fortran
character(:), allocatable :: hostname()
```

Get username of the current user.

```fortran
character(:), allocatable :: get_username()
```

Get owner name of file or directory.

```fortran
character(:), allocatable :: get_owner_name("my/file.txt")
```

Get owner group of file or directory.

```fortran
character(:), allocatable :: get_owner_group("my/file.txt")
```

Get full path of main executable, regardless of current working directory

```fortran
character(:), allocatable :: exe_path()
```

Get full path of **SHARED LIBRARY**, regardless of current working directory.
If static library, returns empty string.
To use `lib_path()`, build Ffilesystem with `cmake -DBUILD_SHARED_LIBS=on`

```fortran
character(:), allocatable :: lib_path()
```

Get current working directory

```fortran
character(:), allocatable :: get_cwd()
```

Change current working directory (chdir):

```fortran
logical :: ok
ok = set_cwd("my/path")
```

Get command argument (allocatable character function vs. Fortran 2003 subroutine get_command_argument()).
If the argument does not exist, returns empty string.

```fortran
character(:), allocatable :: get_arg(i)
```

---

Get environment variable (allocatable character function vs. Fortran 2003 subroutine get_environment variable()):

```fortran
character(:), allocatable :: getenv(name)
```

Set environment variable:

```fortran
logical, optional :: ok

call setenv(name, value, ok=ok)
```

Get system or user temporary directory:

```fortran
character(:), allocatable :: get_tempdir()
```

Make a path absolute if relative:

```fortran
function absolute(path, base)
!! if path is absolute, return path
!! if path is relative, base path
!!
!! idempotent iff top_path is absolute

character(:), allocatable :: absolute
character(*), intent(in) :: path, base
```

Tell characteristics of the computing platform such as operating system:

```fortran
! logical based on C preprocessor

is_admin()
is_bsd()
is_unix()
is_linux()
is_windows()
is_macos()
is_cygwin()
is_mingw()
is_msvc()
```

---

```fortran
logical :: is_rosetta()
```

Tell if program is running via macOS Rosetta 2 translation (macOS Apple Silicon running Intel binaries).

---

logical: is the user running as admin / root / superuser:

```fortran
is_admin()
```

```fortran
C_INT  is_wsl()  !< Windows Subsystem for Linux > 0 if true
```

We have "fs_getpid()" as the function name to avoid clashes with proprietary compiler function.

```fortran
fs_getpid()  !< process ID
```

```fortran
