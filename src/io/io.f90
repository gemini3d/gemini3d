module io
!! HANDLES INPUT AND OUTPUT OF PLASMA STATE PARAMETERS (NOT GRID INPUTS)
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use, intrinsic :: iso_c_binding, only: c_int
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

use phys_consts, only : kB,ms,pi,lsp,wp,lwave
use pathlib, only: mkdir, copyfile
use mpimod, only: bcast_recv, bcast_send, gather_send, gather_recv,  myid, &
  tagns, tagvs1, tagv2, tagv3, tagAur, tagTs, tagJ1, tagJ2, tagJ3
use grid, only : gridflag,flagswap,lx1,lx2,lx3,lx2all, lx3all

implicit none

private
public :: create_outdir, &
  input_plasma, output_plasma, input_plasma_currents, &
  create_outdir_mag, output_magfields, &
  create_outdir_aur, output_aur

interface ! aurora.f90
module subroutine create_outdir_aur(outdir)
character(*), intent(in) :: outdir
end subroutine create_outdir_aur

module subroutine output_aur(outdir,flagglow,ymd,UTsec,iver, out_format)
character(*), intent(in) :: outdir, out_format
integer, intent(in) :: flagglow
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur

module subroutine output_aur_workers(iver)
real(wp), dimension(:,:,:), intent(in) :: iver
end subroutine output_aur_workers
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

module subroutine input_plasma(x1,x2,x3all,indatsize,indatfile,ns,vs1,Ts)
real(wp), dimension(-1:), intent(in) :: x1, x2, x3all
character(*), intent(in) :: indatsize, indatfile
real(wp), dimension(-1:,-1:,-1:,:), intent(out) :: ns,vs1,Ts
end subroutine input_plasma

module subroutine input_plasma_currents(outdir,out_format,flagoutput,ymd,UTsec,J1,J2,J3)
character(*), intent(in) :: outdir, out_format
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(:,:,:), intent(out) :: J1,J2,J3
end subroutine input_plasma_currents

module subroutine output_plasma(outdir,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3, out_format)
character(*), intent(in) :: outdir, out_format
integer, intent(in) :: flagoutput
integer, dimension(3), intent(in) :: ymd
real(wp), intent(in) :: UTsec
real(wp), dimension(-1:,-1:,-1:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
real(wp), dimension(:,:,:), allocatable, intent(inout) :: Phiall     !these jokers may not be allocated, but this is allowed as of f2003
real(wp), dimension(:,:,:), intent(in) :: J1,J2,J3
end subroutine output_plasma

end interface


interface ! output.f90
module subroutine create_outdir(outdir,infile,indatsize,indatgrid,indatfile, &
  flagdneu,sourcedir,flagprecfile,precdir,flagE0file,E0dir)
character(*), intent(in) :: outdir, & !< command line argument output directory
  infile, & !< command line argument input file
  indatsize,indatgrid,indatfile,sourcedir, precdir,E0dir
integer, intent(in) :: flagdneu, flagprecfile, flagE0file
end subroutine create_outdir
end interface


end module io
