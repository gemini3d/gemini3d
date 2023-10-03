#include <iostream>
#include <cstdlib>

#include "gemini3d.h"

void fluid_adv(double*, double*, int*, double*, bool*, int*, int*, double*, double*, double*, int*, void*, void*, void*);

// top-level module calls for gemini simulation
int gemini_main(struct params* ps, int* plid2in, int* plid3in){

  int lx1,lx2,lx3;
  int lx2all,lx3all;
  int lsp;
  double UTsec;
  int ymd[3];
  double* fluidvars;
  double* fluidauxvars;
  double* electrovars;    // pointers modifiable by fortran
  double t=0.0, dt=1e-4;
  double tout, tneuBG, tglowout, tdur, tmilestone=0;
  int it, iupdate;
  int flagoutput;
  bool first,flagneuBG;
  int flagdneu;
  double dtneu,dtneuBG;
  int myid,lid;
  void* cfgC;
  void* intvars;
  void* xC;
  int xtype;

  /* Basic setup */
  mpisetup_C();                               // organize mpi workers
  mpiparms_C(&myid,&lid);                     // information about worker number, etc.

  /* Command line and config structure setup */
  // cli_config_gridsize_C(ps, plid2in, plid3in, &cfgC);    // handling of input data, create internal fortran type with parameters for run

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Allocations happen in this block
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  gemini_cfg_alloc_C(&cfgC);
  cli_in_C(ps,plid2in,plid3in,&cfgC);
  read_config_in_C(ps,&cfgC);
  grid_size_in_C(&cfgC);

  // Grab some variables out of fortran modules
  get_fullgrid_size_C(&lx1,&lx2all,&lx3all);                        // retrieve sizes that are stored in the grid module
  init_procgrid_C(&lx2all,&lx3all,plid2in,plid3in);                 // compute process grid for this run
  get_config_vars_C(&cfgC, &flagneuBG,&flagdneu,&dtneuBG,&dtneu);   // export config type properties as C variables, for use in main

  // Once the process grid is set we can compute subgrid sizes
  calc_subgrid_size_in_C(&lx2all,&lx3all);

  /* Main needs to know the grid sizes and species numbers */
  get_subgrid_size_C(&lx1,&lx2,&lx3);     // once grid is input we need to know the subgrid sizes based on no of workers and overall size
  get_species_size_C(&lsp);               // so main knows the number of species used

  // Allocate memory and get pointers to blocks of data
  //gemini_alloc(&fluidvars,&fluidauxvars,&electrovars);    // allocate space in fortran modules for data
  std::cout << "start C allocations:  " << lx1 << " " << lx2 << " " << lx3 << std::endl;
  fluidvars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*5*lsp*sizeof(double));
  fluidauxvars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*(2*lsp+9)*sizeof(double));
  electrovars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*7*sizeof(double));
  if (! fluidvars){
    std::cerr << "fluidvars failed malloc, worker: " << myid << std::endl;
    return 1;
  }
  if (! fluidauxvars){
    std::cerr << "fluiduxvars failed malloc, worker: " << myid << std::endl;
    return 1;
  }
  if (! electrovars){
    std::cerr << "electrovars failed malloc, worker: " << myid << std::endl;
    return 1;
  }
  gemini_work_alloc_C(&cfgC,&intvars);
  outdir_fullgridvaralloc_C(&cfgC,&intvars,&lx1,&lx2all,&lx3all);          // create output directory and allocate some module space for potential
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /* Get input grid from file */
  read_grid_C(&cfgC, &xtype, &xC);                              // read the input grid from file, storage as fortran module object

  /* initialize state variables from input file */
  get_initial_state_C(&cfgC,&fluidvars,&electrovars,&intvars,&xtype,&xC,&UTsec,&ymd[0],&tdur,&t,&tmilestone);
  set_start_values_auxtimevars_C(&it,&t,&tout,&tglowout,&tneuBG);
  set_start_values_auxvars_C(&xtype,&xC,&fluidauxvars);

  /* initialize other file input data */
  init_Efieldinput_C(&cfgC,&xtype,&xC,&dt,&t,&intvars,&ymd[0],&UTsec);
  pot2perpfield_C(&xtype,&xC,&electrovars);

  BGfield_Lagrangian_C(&cfgC, &xtype, &xC, &electrovars, &intvars);
  init_precipinput_C(&cfgC,&xtype,&xC,&dt,&t,&ymd[0],&UTsec,&intvars);
  msisinit_C(&cfgC);
  init_neutralBG_C(&cfgC,&xtype,&xC,&dt,&t,&ymd[0],&UTsec,&intvars);
  init_neutralperturb_C(&dt,&cfgC,&xtype,&xC,&intvars,&ymd[0],&UTsec);

  /* Compute initial drift velocity */
  get_initial_drifts_C(&cfgC, &xtype, &xC, &fluidvars, &fluidauxvars, &electrovars, &intvars);

  /* Control console printing for, actually superfluous FIXME */
  set_update_cadence_C(&iupdate);

  while(t<tdur){
    dt_select_C(&cfgC,&xtype,&xC,&fluidvars,&fluidauxvars,&it,&t,&tout,&tglowout,&dt);
    if (myid ==0){
      std::cout << " ...Selected time step (seconds) " << dt << std::endl;
    }

    // neutral data
    if (it!=1 && flagneuBG && t>tneuBG){
      neutral_atmos_winds_C(&cfgC,&xtype,&xC,&ymd[0],&UTsec,&intvars);
      neutral_atmos_wind_update_C(&intvars);
      tneuBG+=dtneuBG;
      if (myid==0){
        std::cout << " Computed neutral background..." << std::endl;
      }
    }
    if (flagdneu==1){
      neutral_perturb_C(&cfgC,&intvars,&xtype,&xC,&dt,&t,&ymd[0],&UTsec);
      if (myid==0){
        std::cout << " Computed neutral perturbations..." << std::endl;
      }
    }

    // call electrodynamics solution
    //std::cout << " Start electro solution..." << std::endl;
    electrodynamics_C(&cfgC,&fluidvars,&fluidauxvars,&electrovars,&intvars,&xtype,&xC,&it,&t,&dt,&ymd[0],&UTsec);
    //std::cout << " Computed electrodynamics solutions..." << std::endl;

    // advance the fluid state variables
    first=it==1;
    fluid_adv(&t,&dt,&ymd[0],&UTsec,&first,&lsp,&myid,fluidvars,fluidauxvars,electrovars,&xtype,cfgC,xC, intvars);
    //std::cout << " Computed fluid update..." << std::endl;

    check_finite_output_C(&cfgC,&fluidvars,&electrovars,&t);
    it+=1; t+=dt;
    dateinc_C(&dt,&ymd[0],&UTsec);
    check_dryrun_C(&cfgC);
    check_fileoutput_C(&cfgC,&fluidvars,&electrovars,&intvars,&t,&tout,&tglowout,&tmilestone,&flagoutput,&ymd[0],&UTsec);
    if (myid==0){
      std::cout << " Time step " << it << " finished: " << ymd[0] << " " << ymd[1] << " " << ymd[2] << " " << UTsec << " " << t << std::endl;
      //std::cout << " Output cadence variables:  " << tout << " " << tglowout << " " << tmilestone << std::endl;
    }
  }

  /* Call deallocation procedures */
  clear_dneu_C(&intvars);
  free(fluidvars); free(fluidauxvars); free(electrovars);
  gemini_work_dealloc_C(&cfgC,&intvars);
  gemini_cfg_dealloc_C(&cfgC);
  return 0;
}


/* note that none of the pointer locations will be modified, e.g. with malloc, etc. so these pointers can be passed by value */
void fluid_adv(double* pt, double* pdt, int* pymd, double* pUTsec, bool* pfirst, int* plsp, int* pmyid,
  double* fluidvars, double* fluidauxvars, double* electrovars,
  int* pxtype,
  void* cfgC, void* xC, void* intvars)
  {
  double f107,f107a;
  double gavg,Tninf;

  int one=1,two=2,three=3;    // silly but I need some way to pass these ints by reference to fortran...

  /* Set up variables for the time step */
  get_solar_indices_C(&cfgC, &f107,&f107a);   // FIXME: do we really need to return the indices???
  v12rhov1_C(&fluidvars,&fluidauxvars);
  T2rhoe_C(&fluidvars,&fluidauxvars);

  /* Advection substep */
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This old haloing code probably has no real benefit except for not haloing both cells of hte velocities
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // halo_interface_vels_allspec_C(pxtype,&xC,&fluidvars,plsp);
  // interface_vels_allspec_C(&fluidvars,&intvars,plsp);
  // set_global_boundaries_allspec_C(pxtype,&xC,&fluidvars,&fluidauxvars,&intvars,plsp);
  // halo_allparams_C(pxtype, &xC, &fluidvars, &fluidauxvars);

  // Probably very little drawback to doing things this more general way
  set_global_boundaries_allspec_C(pxtype,&xC,&fluidvars,&fluidauxvars,&intvars,plsp);
  halo_fluidvars_C(pxtype, &xC, &fluidvars, &fluidauxvars);
  interface_vels_allspec_C(&fluidvars,&intvars,plsp);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sweep3_allparams_C(&fluidvars,&fluidauxvars,&intvars,pxtype,&xC,pdt);
  sweep1_allparams_C(&fluidvars,&fluidauxvars,&intvars,pxtype,&xC,pdt);
  halo_allparams_C(pxtype, &xC, &fluidvars, &fluidauxvars);
  sweep2_allparams_C(&fluidvars,&fluidauxvars,&intvars,pxtype,&xC,pdt);
  rhov12v1_C(&fluidvars,&fluidauxvars);
  clean_param_C(&one, pxtype, &xC, &fluidvars);
  clean_param_C(&two, pxtype, &xC, &fluidvars);

  /* Compression substep */
  VNRicht_artvisc_C(&fluidvars,&intvars);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // This old haloing code does have benefit since it doesn't automatically halo everything.  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RK2_prep_mpi_allspec_C(pxtype, &xC, &fluidvars);

  // Substantial drawback here because we are unneccessary haloing
  halo_fluidvars_C(pxtype,&xC,&fluidvars,&fluidauxvars);
  RK2_global_boundary_allspec_C(pxtype, &xC, &fluidvars);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  compression_C(&fluidvars,&fluidauxvars,&intvars,pxtype,&xC,pdt);
  rhoe2T_C(&fluidvars,&fluidauxvars);
  clean_param_C(&three, pxtype, &xC, &fluidvars);

  /* Energy diffusion substep */
  energy_diffusion_C(&cfgC,pxtype,&xC,&fluidvars,&electrovars,&intvars,pdt);
  clean_param_C(&three, pxtype, &xC, &fluidvars);
  T2rhoe_C(&fluidvars,&fluidauxvars);

  /* Prep for sources step - all workers must have a common average gravity and exospheric temperature */
  get_gavg_Tinf_C(&intvars, &gavg ,&Tninf);
  /* Sources substep and finalize solution for this time step */
  source_loss_allparams_C(&cfgC,&fluidvars,&fluidauxvars,&electrovars,&intvars,pxtype,&xC,pdt,pt,pymd,pUTsec,&f107a,&f107,pfirst,&gavg,&Tninf);    // note that this includes and conversion of internal energy density and momentum density back to temp and veloc...
  clean_param_C(&three, pxtype, &xC, &fluidvars);
  clean_param_C(&two, pxtype, &xC, &fluidvars);
  clean_param_C(&one, pxtype, &xC, &fluidvars);

  // Fix electron veloc???
}
