module gemini_init

use gemini3d_config, only : gemini_cfg
use mpimod, only : mpi_cfg
use filesystem, only : assert_is_file, assert_is_dir
use phys_consts, only: mindens, mindensnull, mindensdiv

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

  ! base
  print '(A)', 'simulation directory: ' // cfg%outdir
  print '(A51,I6,A1,I0.2,A1,I0.2)', ' start year-month-day:  ', cfg%ymd0(1), '-', cfg%ymd0(2),'-', cfg%ymd0(3)
  print '(A51,F10.3)', 'start time:  ',cfg%UTsec0
  print '(A51,F10.3)', 'duration:  ',cfg%tdur
  print '(A51,F10.3)', 'output every:  ',cfg%dtout
  print*, 'F10.7 and geomagnetic indices:  ',cfg%activ
  print*, 'Top boundary electron temperature:  ',cfg%Teinf 

  ! flags
  print*, 'Potential solve:  ', cfg%potsolve
  print*, 'Periodic:  ', cfg%flagperiodic
  print*, 'Output type:  ', cfg%flagoutput

  ! files
  print '(A,/,A,/,A,/,A)', 'gemini.f90: using input data files:', cfg%indatsize, cfg%indatgrid, cfg%indatfile

  ! neutral perturb
  if(cfg%flagdneu==1) then
    if (.not. (cfg%interptype==5 .or. cfg%interptype==6)) call assert_is_dir(cfg%sourcedir)
    print *, 'Neutral disturbance mlat,mlon:  ',cfg%sourcemlat,cfg%sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',cfg%dtneu
    print *, 'Neutral grid resolution (m):  ',cfg%drhon,cfg%dzn
    print *, 'Neutral disturbance data files located in directory:  ',cfg%sourcedir
  else
    print *, "no neutral disturbance specified."
  end if

  ! precip
  if (cfg%flagprecfile==1) then
    call assert_is_dir(cfg%precdir)
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',cfg%dtprec
    print *, 'Precipitation file input source directory:  ' // cfg%precdir
  else
    print *, "no precipitation specified"
  end if

  ! efield
  if(cfg%flagE0file==1) then
    call assert_is_dir(cfg%E0dir)
    print *, 'Electric field file input cadence (s):  ',cfg%dtE0
    print *, 'Electric field file input source directory:  ' // cfg%E0dir
  else
    print *, "no Efield specified"
  end if

  ! solflux
  if(cfg%flagsolfluxfile==1) then
    call assert_is_dir(cfg%solfluxdir)
    print *, 'Solar flux file input cadence (s):  ',cfg%dtsolflux
    print *, 'Solar flux file input source directory:  ' // cfg%solfluxdir
  else
    print *, "no solar flux specified"
  end if

  ! neutral_BG
  if (cfg%flagneuBG) then
    print*, 'Neutral background updated at:  ',cfg%dtneuBG
    print*, 'Using MSIS version:  ', cfg%msis_version
  end if

  ! neutral_BG file
  if (cfg%flagneutralBGfile==1) then
    print*, 'Neutral background file input cadence (s):  ', cfg%dtneuBGfile
    print*, 'Neutral background file input directory:  ',cfg%neutralBGdir
  end if

  ! glow
  if (cfg%flagglow==1) then
    print *, 'GLOW enabled for auroral emission calculations.'
    print *, 'GLOW electron transport calculation cadence (s): ', cfg%dtglow
    print *, 'GLOW auroral emission output cadence (s): ', cfg%dtglowout
  else
    print *, "GLOW disabled"
  end if

  ! EIA
  if (cfg%flagEIA) then
    print*, 'EIA enables with peok equatorial drift:  ',cfg%v0equator
  else
    print*, 'EIA disabled'
  end if

  ! precip_BG
  print*, 'Background precipitation has total energy flux and energy:  ',cfg%PhiWBG,cfg%W0BG

  ! Jpar
  if (cfg%flagJpar) then
    print*, 'Parallel current calculation enabled.'
  else
    print*, 'Parallel current calculation disabled.'
  end if

  ! capacitance
  print*, 'Inertial capacitance calculation type:  ',cfg%flagcap

  ! diffusion
  print*, 'Diffusion solve type:  ',cfg%diffsolvetype

  ! milestone
  if (cfg%mcadence > 0) then
    print*, 'Milestone output selected; cadence (every nth output) of:  ',cfg%mcadence
  else
    print*, 'Milestone output disabled.'
  end if

  ! gravdrift
  if (cfg%flaggravdrift) then
    print*, 'Gravitational drift terms enabled.'
  else
    print*, 'Gravitaional drift terms disabled.'
  end if

  ! lagrangian
  if (cfg%flaglagrangian) then
    print*, 'Lagrangian grid enabled.'
  else
    print*, 'Lagrangian grid disabled'
  end if

  ! diamagnetic
  if (cfg%flagdiamagnetic) then
    print*, 'Diamagnetic drift terms enabled.'
  else
    print*, 'Diamagnetic drift terms disabled.'
  end if

  ! twoway_coupled
  if (cfg%flagtwoway) then
    print*, 'Two-way coupling enabled.'
  else
    print*, 'Two-way coupling disabled.'
  end if

  ! nodivJ0
  if (cfg%flagnodivJ0) then
    print*, 'Excluding background current divergence.'
  else
    print*, 'Including background current divergence.'
  end if

  ! FBI
  if (cfg%flagFBI>=1) then
    print*, 'Anomalous electrojet heating enabled.'
  else
    print*, 'Anomalous electrojet heating disabled.'
  end if
  if (cfg%flagFBI==2) then
    print*, 'Nonlinear current enabled.'
  else
    print*, 'Nonlinear current disabled.'
  end if

  ! evibcool
  if (cfg%flagevibcool==1) then
    print*, 'Using updated model for electron inelastic collisions.'
  else
    print*, 'Using legacy model for electron inelastic collisions.'
  end if

  ! magpole 
  if (cfg%flagmagpole) then
    print*, 'Using year-based magnetic pole.'
  else
    print*, 'Using default magnetic pole.'
  end if

  ! J1ve
  if (cfg%flagJ1ve) then
    print*, 'Computing parallel electron drift from parallel current density.'
  else
    print*, 'Computing parallel electron drift from ambipolar assumption.'
  end if

  ! nightQ
  if (cfg%flagnightQ) then
    print*, 'Using nighttime ionization calculation.'
  else
    print*, 'Neglecting nighttime ionization.'
  end if

  ! fang
  print*, 'Fang impact ionization version:  ',cfg%flag_fang

  ! fang_pars
  if (cfg%flag_fang==0) then
    print*, 'diff_num_flux:  ',cfg%diff_num_flux
    print*, 'kappa:  ', cfg%kappa
    print*, 'bimax_frac:  ', cfg%bimax_frac
    print*, 'W0_char:  ', cfg%W0_char
  end if

  ! mindens_user
  print*, 'Density fill values:  ',mindens,mindensnull,mindensdiv
   
  print *,  '**************** end input config ***************'
end if


end subroutine check_input_files

end module gemini_init
