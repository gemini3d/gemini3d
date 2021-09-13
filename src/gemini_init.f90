module gemini_init

implicit none (type, external)

private
public :: find_config

contains

subroutine find_config(cfg)

type(gemini_cfg), intent(inout) :: cfg

logical :: exists
integer :: i
character(*), parameter :: locs(2) = [character(18) :: "/inputs/config.nml", "/config.nml"]
character(:), allocatable :: loc

do i = 1, size(locs)
  loc = trim(cfg%outdir // locs(i))
  inquire(file=loc, exist=exists)
  if (exists) then
    cfg%infile = loc
    return
  endif
end do

error stop 'GEMINI3D: could not find configuration file (config.nml) in ' // cfg%outdir

end subroutine find_config


end module gemini_init
