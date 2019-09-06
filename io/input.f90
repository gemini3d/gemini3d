submodule (io) input
implicit none
contains

module procedure read_configfile
!! READS THE INPUT CONFIGURAITON FILE ANDE ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
character(256) :: buf, indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir
integer :: u, i
real(wp) :: NaN
logical :: is_nml

namelist /base/ ymd, UTsec0, tdur, dtout, activ, tcfl, Teinf, &
  potsolve, flagperiodic, flagoutput, flagcap, indat_size, indat_grid, indat_file, &
  flagdneu, flagprecfile, flagE0file, flagglow
namelist /neutral_perturb/ interptype, sourcemlat, sourcemlon, dxn, drhon, dzn, source_dir
namelist /precip/ dtprec, prec_dir
namelist /efield/ dtE0, E0_dir
namelist /glow/ dtglow, dtglowout

NaN = ieee_value(0._wp, ieee_quiet_nan)

is_nml = infile(len(infile)-3:len(infile)) == '.nml'

!> READ CONFIG.DAT FILE FOR THIS SIMULATION
open(newunit=u, file=infile, status='old', action='read')
if (is_nml) then
  read(u,nml=base)
  indatsize = expanduser(indat_size)
  indatgrid = expanduser(indat_grid)
  indatfile = expanduser(indat_file)
else
  read(u,*) ymd(3),ymd(2),ymd(1)
  read(u,*) UTsec0
  read(u,*) tdur
  read(u,*) dtout
  read(u,*) activ(1),activ(2),activ(3)
  read(u,*) tcfl
  read(u,*) Teinf
  read(u,*) potsolve
  read(u,*) flagperiodic
  read(u,*) flagoutput
  read(u,*) flagcap  ! line 11 config.ini
  read(u,'(a256)') buf  !! format specifier needed, else it reads just one character
  indatsize = expanduser(buf)
  read(u,'(a256)') buf
  indatgrid = expanduser(buf)
  read(u,'(a256)') buf
  indatfile = expanduser(buf)   ! line 14
endif
!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (myid==0) then
  print '(A,I6,A1,I0.2,A1,I0.2)', infile//': simulation ymd is:  ',ymd(1),'/',ymd(2),'/',ymd(3)
  print '(A51,F10.3)', 'start time is:  ',UTsec0
  print '(A51,F10.3)', 'duration is:  ',tdur
  print '(A51,F10.3)', 'output every:  ',dtout
  print *, '...using input data files:  '
  print *, '  ',indatsize
  print *, '  ',indatgrid
  print *, '  ',indatfile
end if

!> NEUTRAL PERTURBATION INPUT INFORMATION
!> defaults
interptype=0
sourcemlat=0
sourcemlon=0
dtneu=0
drhon=0
dzn=0
dxn=0
if(is_nml) then
  if (flagdneu == 1) then
    read(u, nml=neutral_perturb)
    sourcedir = expanduser(source_dir)
  endif
else
  read(u,*) flagdneu  ! line 15
  if( flagdneu==1) then
    read(u,*) interptype
    read(u,*) sourcemlat,sourcemlon
    read(u,*) dtneu
    if (interptype==3) then     !read in extra dxn if 3D
      read(u,*) dxn,drhon,dzn
    else
      read(u,*) drhon,dzn
    end if
    read(u,'(A256)') buf
    sourcedir = expanduser(buf)
  endif
end if

if(flagdneu==1 .and. myid ==0) then
  print *, 'Neutral disturbance mlat,mlon:  ',sourcemlat,sourcemlon
  print *, 'Neutral disturbance cadence (s):  ',dtneu
  print *, 'Neutral grid resolution (m):  ',drhon,dzn
  print *, 'Neutral disturbance data files located in directory:  ',sourcedir
end if

!> PRECIPITATION FILE INPUT INFORMATION
!> defaults
dtprec=0
if(is_nml) then
  if (flagprecfile == 1) then
    read(u, nml=precip)
    precdir = expanduser(prec_dir)
  endif
else
  read(u,*, iostat=i) flagprecfile
  if (flagprecfile==1) then
  !! get the location of the precipitation input files
    read(u,*) dtprec
    read(u,'(A256)') buf
    precdir = expanduser(buf)
  end if
end if

if (flagprecfile==1 .and. myid==0) then
  print '(A,F10.3)', 'Precipitation file input cadence (s):  ',dtprec
  print *, 'Precipitation file input source directory:  '//precdir
end if

!> ELECTRIC FIELD FILE INPUT INFORMATION
!> defaults
dtE0=0
if(is_nml) then
  if (flagE0file == 1) then
    read(u, nml=efield)
    E0dir = expanduser(E0_dir)
  endif
else
  read(u,*, iostat=i) flagE0file
  if (flagE0file==1) then
  !! get the location of the precipitation input files
    read(u,*) dtE0
    read(u,'(a256)') buf
    E0dir = expanduser(buf)
  end if
end if

if(flagE0file==1 .and. myid==0) then
  print *, 'Electric field file input cadence (s):  ',dtE0
  print *, 'Electric field file input source directory:  '//E0dir
end if

!> GLOW ELECTRON TRANSPORT INFORMATION
!> defaults
dtglow=NaN
dtglowout=NaN
if(is_nml) then
  if (flagglow == 1) read(u, nml=glow)
else
  read(u,*, iostat=i) flagglow
  if (flagglow==1) then
    read(u,*, iostat=i) dtglow
    read(u,*, iostat=i) dtglowout
    if (i /= 0) error stop 'did you include GLOW timing config data in the .ini file?'
  end if
end if

if (flagglow==1 .and. myid == 0) then
  print *, 'GLOW enabled for auroral emission calculations.'
  print *, 'GLOW electron transport calculation cadence (s): ', dtglow
  print *, 'GLOW auroral emission output cadence (s): ', dtglowout
end if

close(u)

end procedure read_configfile

end submodule input