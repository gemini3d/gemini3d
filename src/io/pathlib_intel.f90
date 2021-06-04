submodule (pathlib) pathlib_intel

implicit none (type, external)

contains

module procedure directory_exists

inquire(directory=expanduser(path), exist=exists)

end procedure directory_exists

end submodule pathlib_intel
