module gemini_init

use gemini3d_config, only : gemini_cfg
use mpimod, only : mpi_cfg
use filesystem, only : assert_is_file, assert_is_dir

implicit none (type, external)

private
public :: find_config, check_input_files

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


subroutine check_input_files(cfg)

type(gemini_cfg), intent(in) :: cfg

!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (mpi_cfg%myid==0) then
  call assert_is_file(cfg%indatsize)
  call assert_is_file(cfg%indatgrid)
  call assert_is_file(cfg%indatfile)

  print *, '******************** input config ****************'
  print '(A)', 'simulation directory: ' // cfg%outdir
  print '(A51,I6,A1,I0.2,A1,I0.2)', ' start year-month-day:  ', cfg%ymd0(1), '-', cfg%ymd0(2),'-', cfg%ymd0(3)
  print '(A51,F10.3)', 'start time:  ',cfg%UTsec0
  print '(A51,F10.3)', 'duration:  ',cfg%tdur
  print '(A51,F10.3)', 'output every:  ',cfg%dtout
  print '(A,/,A,/,A,/,A)', 'gemini.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

  if(cfg%flagdneu==1) then
    if (.not. (cfg%interptype==5 .or. cfg%interptype==6)) call assert_is_dir(cfg%sourcedir)
    print *, 'Neutral disturbance mlat,mlon:  ',cfg%sourcemlat,cfg%sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',cfg%dtneu
    print *, 'Neutral grid resolution (m):  ',cfg%drhon,cfg%dzn
    print *, 'Neutral disturbance data files located in directory:  ',cfg%sourcedir
  else
    print *, "no neutral disturbance specified."
  end if

  if (cfg%flagprecfile==1) then
    call assert_is_dir(cfg%precdir)
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',cfg%dtprec
    print *, 'Precipitation file input source directory:  ' // cfg%precdir
  else
    print *, "no precipitation specified"
  end if

  if(cfg%flagE0file==1) then
    call assert_is_dir(cfg%E0dir)
    print *, 'Electric field file input cadence (s):  ',cfg%dtE0
    print *, 'Electric field file input source directory:  ' // cfg%E0dir
  else
    print *, "no Efield specified"
  end if

  if (cfg%flagglow==1) then
    print *, 'GLOW enabled for auroral emission calculations.'
    print *, 'GLOW electron transport calculation cadence (s): ', cfg%dtglow
    print *, 'GLOW auroral emission output cadence (s): ', cfg%dtglowout
  else
    print *, "GLOW disabled"
  end if

  if (cfg%msis_version > 0) then
    print '(A,f3.1,A)', 'MSIS ', real(cfg%msis_version)/10, 'enabled for neutral atmosphere calculations.'
  else
    print '(A)', "MSISE00 enabled for neutral atmosphere calculations."
  end if

  if (cfg%flagEIA) then
    print*, 'EIA enables with peok equatorial drift:  ',cfg%v0equator
  else
    print*, 'EIA disabled'
  end if

  if (cfg%flagneuBG) then
    print*, 'Variable background neutral atmosphere enabled at cadence:  ',cfg%dtneuBG
  else
    print*, 'Variable background neutral atmosphere disabled.'
  end if

  print*, 'Background precipitation has total energy flux and energy:  ',cfg%PhiWBG,cfg%W0BG

  if (cfg%flagJpar) then
    print*, 'Parallel current calculation enabled.'
  else
    print*, 'Parallel current calculation disabled.'
  end if

  print*, 'Inertial capacitance calculation type:  ',cfg%flagcap

  print*, 'Diffusion solve type:  ',cfg%diffsolvetype

  if (cfg%mcadence > 0) then
    print*, 'Milestone output selected; cadence (every nth output) of:  ',cfg%mcadence
  else
    print*, 'Milestone output disabled.'
  end if

  if (cfg%flaggravdrift) then
    print*, 'Gravitational drift terms enabled.'
  else
    print*, 'Gravitaional drift terms disabled.'
  end if

  if (cfg%flaglagrangian) then
    print*, 'Lagrangian grid enabled.'
  else
    print*, 'Lagrangian grid disabled'
  end if

  print *,  '**************** end input config ***************'
end if


end subroutine check_input_files

end module gemini_init
