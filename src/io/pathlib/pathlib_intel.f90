submodule (pathlib) pathlib_intel

implicit none (type, external)

contains

module procedure is_dir
inquire(directory=expanduser(path), exist=is_dir)
end procedure is_dir

end submodule pathlib_intel
