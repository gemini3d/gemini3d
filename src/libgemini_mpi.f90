!> mpi-related gemini functionality; all procedures must be bind(C) to allow C/C++ main programs
module gemini3d_mpi

use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_int, c_bool, c_float, c_double
use phys_consts, only: wp,debug
use mpimod, only: mpi_manualgrid, process_grid_auto, mpi_cfg, mpibreakdown
use meshobj, only: curvmesh
use config, only: gemini_cfg
use io, only: output_plasma,output_aur,find_milestone,input_plasma,create_outdir,create_outdir_aur
use potential_comm, only: get_BGEfields,velocities
use grid_mpi, only: grid_drift
use collisions, only: conductivities
use gemini3d, only: cfg,x

implicit none (type, external)
private
public :: init_procgrid, outdir_fullgridvaralloc, get_initial_state, BGfield_Lagrangian, check_dryrun, check_fileoutput,  &
            get_initial_drifts

interface !libgem_mpi_io.f90
  module subroutine outdir_fullgridvaralloc(Phiall,lx1,lx2all,lx3all) bind(C)
    real(kind=c_double), dimension(:,:,:), allocatable, intent(inout) :: Phiall
    integer(kind=c_int), intent(in) :: lx1,lx2all,lx3all
  end subroutine outdir_fullgridvaralloc
  module subroutine get_initial_state(ns,vs1,Ts,Phi,Phiall,UTsec,ymd,tdur) bind(C)
    real(kind=c_double), dimension(:,:,:,:), intent(inout) :: ns,vs1,Ts
    real(kind=c_double), dimension(:,:,:), intent(inout) :: Phi,Phiall
    real(kind=c_double), intent(inout) :: UTsec
    integer(kind=c_int), dimension(3), intent(inout) :: ymd
    real(kind=c_double), intent(inout) :: tdur
  end subroutine get_initial_state
  module subroutine check_fileoutput(t,tout,tglowout,tmilestone,flagoutput,ymd,UTsec,vs2,vs3,ns,vs1,Ts,Phiall,J1,J2,J3,iver) bind(C)
    real(kind=c_double), intent(in) :: t
    real(kind=c_double), intent(inout) :: tout,tglowout,tmilestone
    integer(kind=c_int), intent(inout) :: flagoutput
    integer(kind=c_int), dimension(3), intent(in) :: ymd
    real(kind=c_double), intent(in) :: UTsec
    real(kind=c_double), dimension(:,:,:,:), intent(in) :: vs2,vs3,ns,vs1,Ts
    real(kind=c_double), dimension(:,:,:), allocatable, intent(inout) :: Phiall
    real(kind=c_double), dimension(:,:,:), intent(in) :: J1,J2,J3
    real(kind=c_double), dimension(:,:,:), intent(in) :: iver
  end subroutine
  module subroutine check_dryrun() bind(C)
  end subroutine check_dryrun
end interface

interface !libgem_mpi_drifts.f90
  module subroutine BGfield_Lagrangian(v2grid,v3grid,E1,E2,E3) bind(C)
    real(kind=c_double), intent(inout) :: v2grid,v3grid
    real(kind=c_double), dimension(:,:,:), intent(inout) :: E1,E2,E3
  end subroutine BGfield_Lagrangian
  module subroutine get_initial_drifts(nn,Tn,vn1,vn2,vn3,ns,Ts,vs1,vs2,vs3,B1,E2,E3) bind(C)
    real(kind=c_double), dimension(:,:,:,:), intent(in) :: nn
    real(kind=c_double), dimension(:,:,:), intent(in) :: Tn,vn1,vn2,vn3
    real(kind=c_double), dimension(:,:,:,:), intent(in) :: ns,Ts,vs1
    real(kind=c_double), dimension(:,:,:,:), intent(inout) :: vs2,vs3
    real(kind=c_double), dimension(:,:,:), intent(in) :: B1
    real(kind=c_double), dimension(:,:,:), intent(in) :: E2,E3
  end subroutine get_initial_drifts
end interface
interface !libgem_mpi_par.f90
  module subroutine init_procgrid(lx2all,lx3all,lid2in,lid3in) bind(C)
    integer(kind=c_int), intent(in) :: lx2all,lx3all,lid2in,lid3in
  end subroutine
end interface

contains

end module gemini3d_mpi
