submodule (io) input

implicit none

interface ! path_exists*.f90
module subroutine assert_directory_exists(path)
character(*), intent(in) :: path
end subroutine assert_directory_exists
end interface

contains

module procedure read_configfile
!! READS THE INPUT CONFIGURAITON FILE, ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.
character(256) :: buf, indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir
integer :: i
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

!> READ CONFIG FILE FOR THIS SIMULATION
rawconfig : block
integer :: u
open(newunit=u, file=infile, status='old', action='read')
if (is_nml) then
  read(u,nml=base, iostat=i)
  call check_nml_io(i, infile)
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

call assert_file_exists(indatsize)
call assert_file_exists(indatgrid)
call assert_file_exists(indatfile)

!> PRINT SOME DIAGNOSIC INFO FROM ROOT
if (myid==0) then
  print '(A,I6,A1,I0.2,A1,I0.2)', infile // ' simulation year-month-day is:  ',ymd(1),'-',ymd(2),'-',ymd(3)
  print '(A51,F10.3)', 'start time is:  ',UTsec0
  print '(A51,F10.3)', 'duration is:  ',tdur
  print '(A51,F10.3)', 'output every:  ',dtout
  print '(A,/,A,/,A,/,A)', 'using input data files:', indatsize, indatgrid, indatfile
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
    read(u, nml=neutral_perturb, iostat=i)
    call check_nml_io(i, infile)
    sourcedir = expanduser(source_dir)
    call assert_directory_exists(sourcedir)
  endif
else
  read(u,*, iostat=i) flagdneu  ! line 15
  call check_ini_io(i, infile)
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
    call assert_directory_exists(sourcedir)
  endif
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagdneu/=1) sourcedir = ""

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
    read(u, nml=precip, iostat=i)
    call check_nml_io(i, infile)
    precdir = expanduser(prec_dir)
    call assert_directory_exists(precdir)
  endif
else
  read(u,*, iostat=i) flagprecfile
  call check_ini_io(i, infile)
  if (flagprecfile==1) then
  !! get the location of the precipitation input files
    read(u,*, iostat=i) dtprec
    read(u,'(A256)', iostat=i) buf
    call check_nml_io(i, infile)
    precdir = expanduser(buf)
    call assert_directory_exists(precdir)
  end if
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagprecfile/=1) precdir = ""

if (flagprecfile==1 .and. myid==0) then
  print '(A,F10.3)', 'Precipitation file input cadence (s):  ',dtprec
  print *, 'Precipitation file input source directory:  ' // precdir
end if

!> ELECTRIC FIELD FILE INPUT INFORMATION
!> defaults
dtE0=0
if(is_nml) then
  if (flagE0file == 1) then
    read(u, nml=efield, iostat=i)
    call check_nml_io(i, infile)
    E0dir = expanduser(E0_dir)
    call assert_directory_exists(E0dir)
  endif
else
  read(u,*, iostat=i) flagE0file
  call check_ini_io(i, infile)
  if (flagE0file==1) then
  !! get the location of the precipitation input files
    read(u,*, iostat=i) dtE0
    read(u,'(a256)', iostat=i) buf
    call check_nml_io(i, infile)
    E0dir = expanduser(buf)
    call assert_directory_exists(E0dir)
  end if
end if
!> have to allocate, even when not used, to avoid runtime errors with pickier compilers
if (flagE0file/=1) E0dir = ""

if(flagE0file==1 .and. myid==0) then
  print *, 'Electric field file input cadence (s):  ',dtE0
  print *, 'Electric field file input source directory:  ' // E0dir
end if

!> GLOW ELECTRON TRANSPORT INFORMATION
!> defaults
dtglow=NaN
dtglowout=NaN
if(is_nml) then
  if (flagglow == 1) then
    read(u, nml=glow, iostat=i)
    call check_nml_io(i, infile)
  endif
else
  read(u,*, iostat=i) flagglow
  call check_ini_io(i, infile)
  if (flagglow==1) then
    read(u,*, iostat=i) dtglow
    read(u,*, iostat=i) dtglowout
    call check_nml_io(i, infile)
  end if
end if

if (flagglow==1 .and. myid == 0) then
  print *, 'GLOW enabled for auroral emission calculations.'
  print *, 'GLOW electron transport calculation cadence (s): ', dtglow
  print *, 'GLOW auroral emission output cadence (s): ', dtglowout
end if

close(u)
end block rawconfig

end procedure read_configfile


subroutine check_nml_io(i, filename)
!! checks for EOF and gives helpful error
!! this accomodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename

if (is_iostat_end(i)) then
  write(stderr,*) 'ERROR: ensure there is a trailing blank line in ' // filename
  error stop
endif

if (i /= 0) then
  write(stderr,*) 'ERROR: problem reading ' // filename
  error stop
endif

end subroutine check_nml_io


subroutine check_ini_io(i, filename)
!! checks for EOF and gives helpful error
!! this accomodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename

if (is_iostat_end(i)) return

if (i /= 0) then
write(stderr,*) 'ERROR: problem reading ' // filename
error stop
endif

end subroutine check_ini_io


subroutine assert_file_exists(path)
!! throw error if file does not exist
!! this accomodates non-Fortran 2018 error stop with variable character

character(*), intent(in) :: path
logical :: exists

inquire(file=path, exist=exists)

if (.not.exists) then
  write(stderr,*) 'ERROR: file does not exist' // path
  error stop
endif

end subroutine assert_file_exists

end submodule input
