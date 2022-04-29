!> This module contains C/CXX wrappers for functions in libgemini_mpi.  These routines match those in libgemini_mpi.f90 and are
!!   principally meant to convert the C pointers to various data objects into fortran pointers (including in the case of the
!!   grid a class pointer (pointer to polymorphic object).  Other polymorhpic objects (neutraldata, etc.) are kept in a static
!!   derived type (intvars::gemini_work) and don't need to be passes as class pointers.
module gemini3d_mpi_C

use, intrinsic :: iso_c_binding, only : c_f_pointer, c_ptr, C_INT

use meshobj, only: curvmesh
use config, only: gemini_cfg
use io, only: output_plasma,output_aur,find_milestone,input_plasma,create_outdir
use potential_comm, only: get_BGEfields,velocities
use grid, only: lx1,lx2,lx3
use grid_mpi, only: grid_drift, read_grid
use collisions, only: conductivities
use potentialBCs_mumps, only: init_Efieldinput
use potential_comm,only : pot2perpfield, electrodynamics
use neutral_perturbations, only: init_neutralperturb,neutral_denstemp_update,neutral_wind_update,neutral_perturb
use temporal, only : dt_comm
use sanity_check, only : check_finite_pertub, check_finite_output
use advec_mpi, only: set_global_boundaries_allspec, halo_interface_vels_allspec
use multifluid_mpi, only: halo_allparams
use sources_mpi, only: RK2_prep_mpi_allspec
use ionization_mpi, only: get_gavg_Tinf
use neutral_perturbations, only: clear_dneu

use gemini3d, only: fluidvar_pointers,fluidauxvar_pointers, electrovar_pointers, gemini_work
use gemini3d_mpi, only: outdir_fullgridvaralloc

implicit none (type, external)

contains

  !> create output directory and allocate full grid potential storage
  subroutine outdir_fullgridvaralloc_C(cfgC,intvarsC,lx1,lx2all,lx3all) bind(C)
    type(c_ptr), intent(in) :: cfgC
    type(c_ptr), intent(inout) :: intvarsC
    integer(C_INT), intent(in) :: lx1,lx2all,lx3all

    type(gemini_cfg), pointer :: cfg
    type(gemini_work), pointer :: intvars

    call c_f_pointer(cfgC, cfg)
    call c_f_pointer(intvarsC,intvars)

    call outdir_fullgridvaralloc(cfg, intvars, lx1, lx2all, lx3all)
  end subroutine outdir_fullgridvaralloc_C

end module gemini3d_mpi_C
