submodule(config) config_nml

implicit none (type, external)

contains


module procedure read_nml
!! Reads simulation configuration file in .nml
!! Note that it is best to rewind the file before any read operation, otherwise if the file pointer is already
!! past the group of interest it will (may?) miss that group and return junk.

integer :: u, i
logical :: exists

integer :: ymd(3)
real(wp) :: UTsec0
real(wp) :: tdur
real(wp) :: dtout
real(wp) :: activ(3)
real(wp) :: tcfl
real(wp) :: Teinf
integer :: potsolve, flagperiodic=0, flagoutput, flagcap=0
integer :: interptype
real(wp) :: sourcemlat,sourcemlon
real(wp) :: dtneu
real(wp) :: dxn,drhon,dzn
real(wp) :: dtprec=0
character(256) :: indat_size, indat_grid, indat_file, source_dir, prec_dir, E0_dir
character(4) :: file_format=""  !< need to initialize blank or random invisible fouls len_trim>0
real(wp) :: dtE0=0
integer :: flagdneu, flagprecfile, flagE0file, flagglow !< FIXME: these parameters are ignored, kept temporarily
real(wp) :: dtglow=0, dtglowout=0
logical :: flagEIA
real(wp) :: v0equator
logical :: flagneuBG
real(wp) :: dtneuBG
real(wp) :: PhiWBG,W0BG
logical :: flagJpar
logical :: flgcap
real(wp) :: magcap
integer :: diffsolvetype
integer :: mcadence

namelist /base/ ymd, UTsec0, tdur, dtout, activ, tcfl, Teinf
namelist /files/ file_format, indat_size, indat_grid, indat_file
namelist /flags/ potsolve, flagperiodic, flagoutput, &
flagcap,flagdneu, flagprecfile, flagE0file, flagglow !< FIXME: these last parameters are ignored, kept temporarily for compatibility, should be removed
namelist /neutral_perturb/ flagdneu, interptype, sourcemlat, sourcemlon, dtneu, dxn, drhon, dzn, source_dir
namelist /precip/ dtprec, prec_dir
namelist /efield/ dtE0, E0_dir
namelist /glow/ dtglow, dtglowout
namelist /EIA/ flagEIA,v0equator
namelist /neutral_BG/ flagneuBG,dtneuBG
namelist /precip_BG/ PhiWBG,W0BG
namelist /Jpar/ flagJpar
namelist /capacitance/ flagcap,magcap     ! later need to regroup these in a way that is more logical now there are so many more inputs
namelist /diffusion/ diffsolvetype
namelist /milestone/ mcadence

open(newunit=u, file=cfg%infile, status='old', action='read')

read(u, nml=base, iostat=i)
call check_nml_io(i, cfg%infile, "base")
cfg%ymd0 = ymd
cfg%UTsec0 = UTsec0
cfg%tdur = tdur
cfg%dtout = dtout
cfg%activ = activ
cfg%tcfl = tcfl
cfg%Teinf = Teinf

rewind(u)
read(u, nml=flags, iostat=i)
call check_nml_io(i, cfg%infile, "flags")
cfg%potsolve = potsolve
cfg%flagperiodic = flagperiodic
cfg%flagoutput = flagoutput

rewind(u)
read(u, nml=files, iostat=i)
call check_nml_io(i, cfg%infile, "files")

!> auto file_format if not specified
if (len_trim(file_format) > 0) then
  cfg%out_format = trim(file_format)
else
  file_format = get_suffix(indat_size)
  cfg%out_format = file_format(2:)
endif

!> absolute paths or paths relative to cfg%outdir
cfg%indatsize = make_absolute(indat_size, cfg%outdir)
cfg%indatgrid = make_absolute(indat_grid, cfg%outdir)
cfg%indatfile = make_absolute(indat_file, cfg%outdir)

if (namelist_exists(u, "neutral_perturb", verbose)) then
  cfg%flagdneu = 1
  rewind(u)
  read(u, nml=neutral_perturb, iostat=i)
  call check_nml_io(i, cfg%infile, "neutral_perturb")
  cfg%sourcedir = make_absolute(source_dir, cfg%outdir)
  cfg%interptype = interptype
  cfg%sourcemlat = sourcemlat
  cfg%sourcemlon = sourcemlon
  cfg%dtneu = dtneu
  cfg%drhon = drhon
  cfg%dzn = dzn
  cfg%dxn = dxn
else
  cfg%flagdneu = 0
  cfg%sourcedir = ""
endif

if (namelist_exists(u, "precip", verbose)) then
  cfg%flagprecfile = 1
  rewind(u)
  read(u, nml=precip, iostat=i)
  call check_nml_io(i, cfg%infile, "precip")
  cfg%precdir = make_absolute(prec_dir, cfg%outdir)
  cfg%dtprec = dtprec
else
  cfg%flagprecfile = 0
  cfg%precdir = ""
endif

if (namelist_exists(u, "efield", verbose)) then
  cfg%flagE0file = 1
  rewind(u)
  read(u, nml=efield, iostat=i)
  call check_nml_io(i, cfg%infile, "efield")
  cfg%E0dir = make_absolute(E0_dir, cfg%outdir)
  cfg%dtE0 = dtE0
else
  cfg%flagE0file = 0
  cfg%E0dir = ""
endif

if (namelist_exists(u, "glow", verbose)) then
  cfg%flagglow = 1
  rewind(u)
  read(u, nml=glow, iostat=i)
  call check_nml_io(i, cfg%infile, "glow")
  cfg%dtglow = dtglow
  cfg%dtglowout = dtglowout
else
  cfg%flagglow = 0
endif

!> EIA (optional)
if (namelist_exists(u,'EIA')) then
  rewind(u)
  read(u, nml=EIA, iostat=i)
  call check_nml_io(i, cfg%infile, "EIA")
  cfg%flagEIA=flagEIA
  cfg%v0equator=v0equator
end if

!> neural background (optional)
if (namelist_exists(u,'neutral_BG')) then
  rewind(u)
  read(u, nml=neutral_BG, iostat=i)
  call check_nml_io(i, cfg%infile, "neutral_BG")
  cfg%flagneuBG=flagneuBG
  cfg%dtneuBG=dtneuBG
end if

!> precip background (optional)
if (namelist_exists(u,'precip_BG')) then
  rewind(u)
  read(u, nml=precip_BG, iostat=i)
  call check_nml_io(i, cfg%infile, "precip_BG")
  cfg%PhiWBG=PhiWBG
  cfg%W0BG=W0BG
end if

!> parallel current density (optional)
if (namelist_exists(u,'Jpar')) then
  rewind(u)
  read(u, nml=Jpar, iostat=i)
  call check_nml_io(i, cfg%infile, "Jpar")
  cfg%flagJpar=flagJpar
end if

!> inertial capacitance (optional)
if (namelist_exists(u,'capacitance')) then
  rewind(u)
  read(u, nml=capacitance, iostat=i)
  call check_nml_io(i, cfg%infile, "capacitance")
  cfg%flagcap=flagcap
  cfg%magcap=magcap
end if

!> diffusion solve type (optional)
if (namelist_exists(u,'diffusion')) then
  rewind(u)
  read(u, nml=diffusion, iostat=i)
  call check_nml_io(i, cfg%infile, "diffusion")
  cfg%diffsolvetype=diffsolvetype
end if

!> information about milestone outputs
if (namelist_exists(u,'milestone')) then
  rewind(u)
  read(u,nml=milestone,iostat=i)
  call check_nml_io(i,cfg%infile,"milestone")
  cfg%mcadence=mcadence
end if

close(u)

end procedure read_nml


logical function namelist_exists(u, nml, verbose)
!! determines if Namelist exists in file

character(*), intent(in) :: nml    ! FIXME:  is it bad to use a keyword as a variable name?
integer, intent(in) :: u
logical, intent(in), optional :: verbose

logical :: debug
integer :: i
character(256) :: line  !< arbitrary length

debug = .false.
if(present(verbose)) debug = verbose

namelist_exists = .false.

rewind(u)

do
  read(u, '(A)', iostat=i) line
  if(i/=0) exit
  if (line(1:1) /= '&') cycle
  if (line(2:) == nml) then
    namelist_exists = .true.
    exit
  end if
end do
rewind(u)

if (debug) print *, 'namelist ', nml, namelist_exists

end function namelist_exists


subroutine check_nml_io(i, filename, namelist)
!! checks for EOF and gives helpful error
!! this accomodates non-Fortran 2018 error stop with variable character

integer, intent(in) :: i
character(*), intent(in) :: filename
character(*), intent(in), optional :: namelist

character(:), allocatable :: nml, msg

if(i==0) return

nml = ""
if(present(namelist)) nml = namelist

if (is_iostat_end(i)) then
  write(stderr,*) 'ERROR: namelist ' // nml // ': ensure there is a trailing blank line in ' // filename
  error stop 5
endif

msg = ""
select case (get_compiler_vendor())
case ("Intel")
  !! https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-list-of-run-time-error-messages
  select case (i)
  case (19)
    msg = "mismatch between variable names in namelist and Fortran code, or problem in variable specification in file"
  case (623)
    msg = "variable specified in Fortran code missing from Namelist file"
  case (17,18,624,625,626,627,628,680,750,759)
    msg = "namelist file format problem"
  end select
case ("GCC", "GNU")
  select case (i)
  case (5010)
    msg = "mismatch between variable names in namelist and Fortran code, or problem in variable specification in file"
  end select
end select


if (len(msg)==0) write(stderr,*) "namelist read error code",i

write(stderr,'(A,/,A)') 'ERROR: reading namelist ', nml, " from ", filename, " problem: ", msg
error stop 5

end subroutine check_nml_io


end submodule config_nml
