submodule (io) input
implicit none
contains

module procedure read_configfile
!! READS THE INPUT CONFIGURAITON FILE ANDE ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
character(256) :: buf
integer :: u, i
real(wp) :: NaN

NaN = ieee_value(0._wp, ieee_quiet_nan)

!> READ CONFIG.DAT FILE FOR THIS SIMULATION
open(newunit=u, file=infile, status='old', action='read')
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
read(u,*) flagcap                 ! line 11 config.ini
read(u,'(a256)') buf
!! format specifier needed, else it reads just one character
indatsize = expanduser(buf)
read(u,'(a256)') buf
indatgrid = expanduser(buf)
read(u,'(a256)') buf
indatfile = expanduser(buf)       ! line 14 config.ini

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
read(u,*, iostat=i) flagdneu
if (i /= 0) then
  write(stderr, *) 'not using perturbation input files since not specified in ',infile
  flagdneu=0
endif
if( flagdneu==1) then
  read(u,*) interptype
  read(u,*) sourcemlat,sourcemlon
  read(u,*) dtneu
  if (interptype==3) then     !read in extra dxn if 3D
    read(u,*) dxn,drhon,dzn
  else
    read(u,*) drhon,dzn
    dxn=0
  end if
  read(u,'(A256)') buf
  sourcedir = expanduser(buf)
  if (myid ==0) then
    print *, 'Neutral disturbance mlat,mlon:  ',sourcemlat,sourcemlon
    print *, 'Neutral disturbance cadence (s):  ',dtneu
    print *, 'Neutral grid resolution (m):  ',drhon,dzn
    print *, 'Neutral disturbance data files located in directory:  ',sourcedir
  end if
else
!! just set it to something
  interptype=0
  sourcemlat=0._wp; sourcemlon=0._wp;
  dtneu=0._wp
  drhon=0._wp; dzn=0._wp;
  sourcedir=''
end if

!> PRECIPITATION FILE INPUT INFORMATION
read(u,*, iostat=i) flagprecfile
if (i /= 0) then
  write(stderr, *) 'not using precipitation input files since not specified in ',infile
  flagprecfile=0
endif
if (flagprecfile==1) then
!! get the location of the precipitation input files
  read(u,*) dtprec

  read(u,'(A256)') buf
  precdir = expanduser(buf)

  if (myid==0) then
    print '(A,F10.3)', 'Precipitation file input cadence (s):  ',dtprec
    print *, 'Precipitation file input source directory:  '//precdir
  end if
else
!! just set results to something
  dtprec=0._wp
  precdir=''
end if

!> ELECTRIC FIELD FILE INPUT INFORMATION
read(u,*, iostat=i) flagE0file
if (i /= 0) then
  write(stderr, *) 'not using E-field input files since not specified in ',infile
  flagE0file=0
endif
if (flagE0file==1) then
!! get the location of the precipitation input files
  read(u,*) dtE0

  read(u,'(a256)') buf
  E0dir = expanduser(buf)

  if (myid==0) then
    print *, 'Electric field file input cadence (s):  ',dtE0
    print *, 'Electric field file input source directory:  '//E0dir
  end if
else                         !just set results to something
  dtE0=0._wp
  E0dir=''
end if

!> GLOW ELECTRON TRANSPORT INFORMATION
read(u,*, iostat=i) flagglow
if (i /= 0) then
  write(stderr, *) 'disabled Glow since not specified in ',infile
  flagglow=0
endif
if (flagglow==1) then
  read(u,*, iostat=i) dtglow
  read(u,*, iostat=i) dtglowout
  if (i /= 0) error stop 'did you include GLOW timing config data in the .ini file?'
  if (myid == 0) then
    print *, 'GLOW enabled for auroral emission calculations.'
    print *, 'GLOW electron transport calculation cadence (s): ', dtglow
    print *, 'GLOW auroral emission output cadence (s): ', dtglowout
  end if
else
  dtglow=NaN
  dtglowout=NaN
end if

close(u)

end procedure read_configfile

end submodule input