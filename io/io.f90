module io
!! HANDLES INPUT AND OUTPUT OF PLASMA STATE PARAMETERS (NOT GRID INPUTS)
use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
use, intrinsic :: iso_c_binding, only: c_int
use phys_consts, only : kB,ms,pi,lsp,wp,lwave
use fsutils, only: expanduser
use std_mkdir, only: mkdir, copyfile
use mpimod, only: bcast_recv, bcast_send, gather_send, gather_recv,  myid, &
  tagns, tagvs1, tagv2, tagv3, tagAur, tagTs, tagJ1, tagJ2, tagJ3
use grid, only : gridflag,flagswap,lx1,lx2,lx3,lx2all, lx3all
use date_formats, only: date_filename
! use logging, only: logger

implicit none

private
public :: read_configfile, create_outdir, &
  input_plasma, output_plasma, input_plasma_currents, &
  create_outdir_mag, output_magfields, &
  create_outdir_aur, output_aur

!> NONE OF THESE VARIABLES SHOULD BE ACCESSED BY PROCEDURES OUTSIDE THIS MODULE
character(:), allocatable, private :: indatfile
!! initial condition data files from input configuration file

interface ! aurora.f90

module subroutine create_outdir_aur(outdir)
character(*), intent(in) :: outdir
end subroutine create_outdir_aur

module subroutine output_aur(outdir,flagglow,ymd,UTsec,iver)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur

module subroutine output_aur_workers(iver)
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_workers

module subroutine output_aur_root(outdir,flagglow,ymd,UTsec,iver)
character(*), intent(in) :: outdir
integer, intent(in) :: flagglow, ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_root

end interface


interface ! mag.f90

module subroutine create_outdir_mag(outdir,fieldpointfile)
character(*), intent(in) :: outdir
character(*), intent(in) :: fieldpointfile
end subroutine create_outdir_mag

module subroutine output_magfields(outdir,ymd,UTsec,Br,Btheta,Bphi)
character(*), intent(in) :: outdir
integer, intent(in) :: ymd(3)
real(wp), intent(in) :: UTsec
real(wp), dimension(:), intent(in)  :: Br,Btheta,Bphi
end subroutine output_magfields

end interface


interface ! plasma.f90

module subroutine input_plasma(x1,x2,x3all,indatsize,ns,vs1,Ts)
real(wp), dimension(-1:), intent(in) :: x1, x2, x3all
character(*), intent(in) :: indatsize
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
end subroutine input_plasma

module subroutine input_plasma_currents(outdir,flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3
end subroutine input_plasma_currents

module subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3)
character(*), intent(in) :: outdir
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall     !these jokers may not be allocated, but this is allowed as of f2003
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
end subroutine output_plasma

end interface

interface ! output.f90
module subroutine create_outdir(outdir,infile,indatsize,indatgrid,flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
character(*), intent(in) :: outdir, & !command line argument output directory
                            infile, & !command line argument input file
                            indatsize,indatgrid,sourcedir, precdir,E0dir
integer, intent(in) :: flagdneu, flagprecfile, flagE0file
end subroutine create_outdir

end interface


contains


subroutine read_configfile(infile,ymd,UTsec0,tdur,dtout,activ,tcfl,Teinf, &
                 potsolve,flagperiodic, flagoutput,flagcap,&
                 indatsize,indatgrid,flagdneu,interptype, &
                 sourcemlat,sourcemlon,dtneu,drhon,dzn,sourcedir,flagprecfile,&
                 dtprec,precdir,flagE0file,dtE0,E0dir,flagglow,dtglow,dtglowout)
!! READS THE INPUT CONFIGURAITON FILE ANDE ASSIGNS VARIABLES FOR FILENAMES, SIZES, ETC.


character(*), intent(in) :: infile
integer, dimension(3), intent(out):: ymd
real(wp), intent(out) :: UTsec0
real(wp), intent(out) :: tdur
real(wp), intent(out) :: dtout
real(wp), dimension(3), intent(out) :: activ
real(wp), intent(out) :: tcfl
real(wp), intent(out) :: Teinf
integer, intent(out) :: potsolve, flagperiodic, flagoutput, flagcap
integer, intent(out) :: flagdneu
integer, intent(out) :: interptype
real(wp), intent(out) :: sourcemlat,sourcemlon
real(wp), intent(out) :: dtneu
real(wp), intent(out) :: drhon,dzn
integer, intent(out) :: flagprecfile
real(wp), intent(out) :: dtprec
character(:), allocatable, intent(out) :: indatsize,indatgrid, sourcedir, precdir, E0dir
integer, intent(out) :: flagE0file
real(wp), intent(out) :: dtE0
integer, intent(out) :: flagglow
real(wp), intent(out) :: dtglow, dtglowout

character(256) :: buf
integer :: u
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
read(u,*) flagcap
read(u,'(a256)') buf
!! format specifier needed, else it reads just one character
indatsize = expanduser(buf)
read(u,'(a256)') buf
indatgrid = expanduser(buf)
read(u,'(a256)') buf
indatfile = expanduser(buf)

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
read(u,*) flagdneu
if( flagdneu==1) then
  read(u,*) interptype
  read(u,*) sourcemlat,sourcemlon
  read(u,*) dtneu
  read(u,*) drhon,dzn
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
read(u,*) flagprecfile
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
read(u,*) flagE0file
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
read(u,*) flagglow
if (flagglow==1) then
  read(u,*) dtglow
  read(u,*) dtglowout
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
end subroutine read_configfile


end module io
