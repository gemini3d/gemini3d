# Boost Filesystem examples

With

* Windows MSVC, GCC 15.1, Clang (ARM or x64)
* Boost::filesystem 1.88.0 or std::filesystem

Observe that Boost::filesystem gives the correct results for \\?\ and \\.\ prefixes, while std::filesystem does not.

```
1: Boost: "C:\a\ball.text"    std::filesystem: "C:\\a\\ball.text"
1:   root_name: "C:"    "C:"
1:   root_path: "C:\"    "C:\\"
1:   root_directory: "\"    "\\"
1: ---------------------
1: Boost: "\\?\C:\a\ball.text"    std::filesystem: "\\\\?\\C:\\a\\ball.text"
1:   root_name: "\\?\C:"    ""
1:   root_path: "\\?\C:\"    "\\"
1:   root_directory: "\"    "\\"
1: ---------------------
1: Boost: "\\.\C:\a\ball.text"    std::filesystem: "\\\\.\\C:\\a\\ball.text"
1:   root_name: "\\.\C:"    ""
1:   root_path: "\\.\C:\"    "\\"
1:   root_directory: "\"    "\\"
1: ---------------------
```
