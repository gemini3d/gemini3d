submodule (pathlib) pathlib_dummy
!! generic routine for non-Intel, non-GCC.
!! beter to make custom per-compiler routine based on pathlib_gcc for other compilers.

implicit none (type, external)

contains

module procedure is_dir
inquire(file=expanduser(path), exist=is_dir)
end procedure is_dir

end submodule pathlib_dummy
