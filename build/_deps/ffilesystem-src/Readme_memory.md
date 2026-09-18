# Dynamic std::string buffers

RAII std::string buffers are used for all C++ string representations of paths, and are automatically freed.
[std::string buffers](https://learn.microsoft.com/en-us/archive/msdn-magazine/2015/july/c-using-stl-strings-at-win32-api-boundaries) are
dynamically sized according to the system maximum path length and the actual path length.
A few functions pre-allocate the maximum system path length using `fs_get_max_path()` because the underlying functions don't return the actual path length.
The length of the path returned is shortened to the actual path length.

These functions are:
