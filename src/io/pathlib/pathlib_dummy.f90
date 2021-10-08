submodule (pathlib) pathlib_dummy
!! generic routine for non-Intel, non-GCC.
!! beter to make custom per-compiler routine based on pathlib_gcc for other compilers.

implicit none (type, external)

contains

module procedure directory_exists

inquire(file=expanduser(path), exist=exists)

end procedure directory_exists

end submodule pathlib_dummy
