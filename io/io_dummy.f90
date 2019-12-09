submodule (pathlib) io_dummy
!! dummy functions

implicit none

contains

module procedure realpath
!! dummy function
!! returns unmodified path

realpath = path

end procedure realpath

end submodule io_dummy