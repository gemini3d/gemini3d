#ifndef GEMINI3D_H
#define GEMINI3D_H

// needed when compiling C files???
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

enum { LMAX = 1000 };

struct params {
  // order and lengths must match in Fortran and C
  // see gemini_main.f90 "cparams"
  bool fortran_nml;
  bool fortran_cli;
  bool debug;
  bool dryrun;
  char out_dir[LMAX];
  // .ini [base]
  int ymd[3];
  float UTsec0;
  float tdur;
  float dtout;
  float activ[3];
  float tcfl;
  float Teinf;
  // .ini
};

// cfgC for now is just passed around as void**, to modify it would have to be struct *


/* This is now housed in the main C program */
extern int gemini_main(struct params *, int*, int*);

extern void help_gemini_bin();

/* interfaces to libgemini */
extern void gemini_cfg_alloc_C(void**);
extern void gemini_cfg_dealloc_C(void**);
extern void cli_in_C(struct params*, int*, int*, void**);
extern void read_config_in_C(struct params*, void**);
extern void grid_size_in_C(void**);
extern void get_fullgrid_size_C(int*, int*, int*);
extern void get_config_vars_C(void**, bool*, int*, double*, double*);
extern void get_subgrid_size_C(int*, int*, int*);
extern void get_species_size_C(int*);
extern void get_fullgrid_lims_C(double*, double*, double*, double*, double*, double*);
extern void gemini_work_alloc_C(void**,void**);
extern void gemini_work_dealloc_C(void**, void**);
//extern void memblock_from_C(double**, double**, double**);
extern void set_start_values_auxtimevars_C(int*, double*, double*, double*, double*);
extern void set_start_timefromcfg_C(void**, int*, double*, double*);
extern void set_start_values_auxvars_C(int*, void**, double**);
extern void get_cfg_timevars_C(void**,double*,bool*,double*,int*,int*);
extern void msisinit_C(void**);
extern void init_neutralBG_C(void**, int*, void**, double*, double*, int*, double*, void**);
extern void set_update_cadence_C(int*);
extern void neutral_atmos_winds_C(void**, int*, void**, int*, double*, void**);
extern void check_finite_output_C(void**, double**, double**, double*);
extern void get_solar_indices_C(void**, double*, double*);
extern void v12rhov1_C(double**, double**);
extern void T2rhoe_C(double**, double**);
extern void interface_vels_allspec_C(double**, void**, int*);
extern void sweep3_allparams_C(double**, double**, void**, int*, void**, double*);
extern void sweep1_allparams_C(double**, double**, void**, int*, void**, double*);
extern void sweep2_allparams_C(double**, double**, void**, int*, void**, double*);
extern void rhov12v1_C(double**, double**);
extern void clean_param_C(int*, int*, void**, double**);
extern void VNRicht_artvisc_C(double**, void**);
extern void compression_C(double**, double**, void**, int*, void**, double*);
extern void rhoe2T_C(double**, double**);
extern void energy_diffusion_C(void**, int*, void**, double**, double**, void**, double*);
extern void source_loss_allparams_C(void**, double**, double**, double**, void**, int*, void**, 
                                      double*, double*, int*, double*, double*, double*, bool*, double*, double*);
extern void check_dryrun_C(void**);
extern void maxcfl_C(double**, int*, void**, double*, double*);
extern void dateinc_C(double*, int*, double*);
extern void interp_file2subgrid_C(void**, int*, void**, double**, double**);
extern void plasma_output_nompi_C(void**, int*, double*, double**, double**, int*, double*, double*, double*);
extern void read_fullsize_gridcenter_C(void**);
extern void grid_from_extents_C(double*,double*,double*,double*,double*,int*,int*,int*,void**);
extern void gemini_grid_alloc_C(double*,double*,double*,int*,int*,int*,int*,void**);
extern void gemini_grid_dealloc_C(int*,void**);
extern void gemini_grid_generate_C(int*, void**);
extern void setv2v3_C(double*, double*);
extern void set_global_boundaries_allspec_C(int*, void**, double**, double**, void**, int*);
extern void checkE1_C(double**, double**, double**, int*);
extern void electrodynamics_test_C(void**,int*,void**,double**,double**,double**,void**);
extern void forceZOH_all_C(double**);
extern void permute_fluidvars_C(double**);
extern void ipermute_fluidvars_C(double**);
extern void tag4refine_C(int*,void**,double**,double**,double**,void**,int*,bool*);
extern void tag4coarsening_C(int*,void**,double**,double**,double**,void**,bool*);
extern void get_grid_magcoords_C(int*,void**,double**,double**,double**);
extern void get_grid_magcoordsi_C(int*,void**,double**,double**,double**);
extern void clean_param_after_regrid_C(int*, int*, void**, double**,void**);
extern void get_locationsi_C(void**,bool*,double*,double*,double*,double**,double**,double**,double**,int*,int*);
extern void get_datainow_ptr_C(void**,double**);
extern void set_datainow_C(void**);
extern void get_neutralperturb_interptype_C(void**,int*);
extern void swap_statevars_C(double**, double**);
extern void interp3_C(double**,double**,double**,int*,int*,int*,double**,double**,int*,int*,
		double*,double*,double*,double**,double**);
extern void interp2_C(double**,double**,int*,int*,double**,double**,int*,int*,double*,double*,double**,
		double**);

/* interfaces for libgemini_mpi */
// some of these will very likely need to be rewritten when used with forestclaw
extern void mpisetup_C();
extern void mpiparms_C(int*, int*);
extern void init_procgrid_C(int*, int*, int*, int*);
extern void read_grid_C(void**, int*, void**);
extern void outdir_fullgridvaralloc_C(void**, void**, int*, int*, int*);
extern void calc_subgrid_size_in_C(int*, int*);
extern void get_initial_state_C(void**, double**, double**, void**, int*, void**, double*, int*, double*, double*, double*);
extern void init_Efieldinput_C(void**, int*, void**, double*, double*, void**, int*, double*);
extern void pot2perpfield_C(int*, void**, double**);
extern void BGfield_Lagrangian_C(void**, int*, void**, double**, void**);
extern void init_precipinput_C(void**, int*, void**, double*, double*, int*, double*, void**);
extern void init_neutralperturb_C(double*, void**, int*, void**, void**, int*, double*);
extern void get_initial_drifts_C(void**, int*, void**, double**, double**, double**, void**);
extern void neutral_atmos_wind_update_C(void**);
extern void neutral_perturb_C(void**, void**, int*, void**, double*, double*, int*, double*);
extern void electrodynamics_C(void**, double**, double**, double**, void**, int*, void**, int*, double*, double*, int*, double*);
extern void get_gavg_Tinf_C(void**, double*, double*);
extern void clear_dneu_C(void**);
extern void halo_interface_vels_allspec_C(int*, void**, double**, int*);
extern void halo_allparams_C(int*, void**, double**, double**);
extern void halo_fluidvars_C(int*, void**, double**, double**);
extern void RK2_prep_mpi_allspec_C(int*, void**, double**);
extern void RK2_global_boundary_allspec_C(int*, void**, double**);
extern void dt_select_C(void**, int*, void**, double**, double**, int*, double*, double*, double*, double*);
extern void check_fileoutput_C(void**, double**, double**, void**, double*, double*, double*, double*, int*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
