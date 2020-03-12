submodule (io:input) path_exists
!! this is for non-Intel compilers
implicit none

contains

module procedure assert_directory_exists
!! throw error if directory does not exist
!! this accomodates non-Fortran 2018 error stop with variable character

logical :: exists

inquire(file=path, exist=exists)

if (.not.exists) then
  write(stderr,*) path // ' directory does not exist'
  error stop
endif

end procedure assert_directory_exists

end submodule path_exists