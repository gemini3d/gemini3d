/* 
This has top-level calls to libgemini in order to run a simulation from C
*/

#include <iostream>
#include <cstdlib>
#include "gemini3d.h"

// top-level module calls for gemini simulation
void gemini_main(struct params* ps, int* plid2in, int* plid3in){
  int ierr;
  int lx1,lx2,lx3;
  int lx2all,lx3all;
  int lsp;
  double UTsec;
  int ymd[3];
  double* fluidvars;
  double* fluidauxvars;
  double* electrovars;    // pointers modifiable by fortran
  double t, dt=1e-6;
  double tout, tneuBG, tglowout, tdur, tmilestone=0;
  double tstart, tfin;
  int it, iupdate;
  int flagoutput;
  double v2grid,v3grid;
  bool first,flagneuBG;
  int flagdneu;
  double dtneu,dtneuBG;
  int myid,lid;
  void fluid_adv(double*, double*, int*, double*, bool*, int*, int*);

  /* Basic setup */
  mpisetup_C();                               // organize mpi workers
  mpiparms_C(&myid,&lid);                     // information about worker number, etc.
  cli_config_gridsize(ps,plid2in,plid3in);    // handling of input data, create internal fortran type with parameters for run
  get_fullgrid_size_C(&lx1,&lx2all,&lx3all);  // read input file that has the grid size information and set it
  init_procgrid(&lx2all,&lx3all,plid2in,plid3in);            // compute process grid for this run
  get_config_vars_C(&flagneuBG,&flagdneu,&dtneuBG,&dtneu);   // export config type properties as C variables, for use in main

  /* Get input grid from file */
  read_grid_C();                              // read the input grid from file, storage as fortran module object

  /* Main needs to know the grid sizes and species numbers */
  get_subgrid_size_C(&lx1,&lx2,&lx3);     // once grid is input we need to know the subgrid sizes based on no of workers and overall size
  get_species_size_C(&lsp);               // so main knows the number of species used

  /* Allocate memory and get pointers to blocks of data */
  //gemini_alloc(&fluidvars,&fluidauxvars,&electrovars);    // allocate space in fortran modules for data
  std::cout << "start C allocations\n";
  fluidvars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*5*lsp*sizeof(double));
  fluidauxvars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*2*lsp*sizeof(double));
  electrovars=(double*) malloc((lx1+4)*(lx2+4)*(lx3+4)*7*sizeof(double));
  if (! fluidvars){
    std::cerr << "fluidvars failed malloc\n";
    return;
  }
  if (! fluidauxvars){
    std::cerr << "fluiduxvars failed malloc\n";
    return;
  }
  if (! electrovars){
    std::cerr << "electrovars failed malloc\n";
    return;
  }
  std::cout << "end C allocations\n";
  memblock_from_C(&fluidvars,&fluidauxvars,&electrovars);
  outdir_fullgridvaralloc(&lx1,&lx2all,&lx3all);          // create output directory and allocate some module space for potential

  /* initialize state variables from input file */
  get_initial_state(&UTsec,&ymd[0],&tdur);
  set_start_values(&it,&t,&tout,&tglowout,&tneuBG);

  /* initialize other file input data */
  std::cout << " Initializing electric field input data..." << std::endl;
  init_Efieldinput_C(&dt,&t,&ymd[0],&UTsec);
  pot2perpfield_C();
  BGfield_Lagrangian(&v2grid,&v3grid);
  std::cout << " Initialize precipitation input data..." << std::endl;
  init_precipinput_C(&dt,&t,&ymd[0],&UTsec);
  std::cout << " Initialize neutral background and input files..." << std::endl;
  msisinit_C();
  init_neutralBG_C(&dt,&t,&ymd[0],&UTsec,&v2grid,&v3grid);
  init_neutralperturb_C(&dt,&ymd[0],&UTsec);

  /* Compute initial drift velocity */
  get_initial_drifts();

  /* Control console printing for, actually superfluous FIXME */
  set_update_cadence(&iupdate);

  while(t<tdur){
    dt_select_C(&it,&t,&tout,&tglowout,&dt);
    if (myid ==0){
      std::cout << " ...Selected time step (seconds) " << dt << std::endl;
    }

    // neutral data
    if (it!=1 && flagneuBG && t>tneuBG){
      neutral_atmos_winds_C(&ymd[0],&UTsec);
      neutral_atmos_wind_update_C(&v2grid,&v3grid);
      tneuBG+=dtneuBG;
      if (myid==0){
        std::cout << " Computed neutral background..." << std::endl;
      }
    }
    if (flagdneu==1){
      neutral_perturb_C(&dt,&t,&ymd[0],&UTsec,&v2grid,&v3grid);
      if (myid==0){
        std::cout << " Computed neutral perturbations..." << std::endl;
      }
    }

    // call electrodynamics solution
    //std::cout << " Start electro solution..." << std::endl;
    electrodynamics_C(&it,&t,&dt,&ymd[0],&UTsec);
    //std::cout << " Computed electrodynamics solutions..." << std::endl;

    // advance the fluid state variables
    first=it==1;
    fluid_adv(&t,&dt,&ymd[0],&UTsec,&first,&lsp,&myid);
    //std::cout << " Computed fluid update..." << std::endl;

    check_finite_output_C(&t);
    it+=1; t+=dt;
    dateinc_C(&dt,&ymd[0],&UTsec);
    check_dryrun();
    check_fileoutput(&t,&tout,&tglowout,&tmilestone,&flagoutput,&ymd[0],&UTsec);
    if (myid==0){
      std::cout << " Time step " << it << " finished: " << ymd[0] << " " << ymd[1] << " " << ymd[2] << " " << UTsec << " " << t << std::endl;
      //std::cout << " Output cadence variables:  " << tout << " " << tglowout << " " << tmilestone << std::endl;
    }
  }

  /* Call deallocation procedures */
  gemini_dealloc(&fluidvars,&fluidauxvars,&electrovars);
  clear_neuBG_C();
  clear_dneu_C();

}

void fluid_adv(double* pt, double* pdt, int* pymd, double* pUTsec, bool* pfirst, int* plsp, int* pmyid){
  double f107,f107a;
  double gavg,Tninf;
  int one=1,two=2,three=3;    // silly but I need some way to pass these ints by reference to fortran...

  /* Set up variables for the time step */
  get_solar_indices_C(&f107,&f107a);   // FIXME: do we really need to return the indices???
  v12rhov1_C();
  T2rhoe_C();

  /* Advection substep */
  halo_interface_vels_allspec_C(plsp);
  interface_vels_allspec_C(plsp);
  set_global_boundaries_allspec_C(plsp);
  halo_allparams_C();
  sweep3_allparams_C(pdt);
  sweep1_allparams_C(pdt);
  halo_allparams_C();
  sweep2_allparams_C(pdt);
  rhov12v1_C();
  clean_param_C(&one);
  clean_param_C(&two);

  /* Compression substep */
  VNRicht_artvisc_C();
  RK2_prep_mpi_allspec_C();
  compression_C(pdt);
  rhoe2T_C();
  clean_param_C(&three);

  /* Energy diffusion substep */
  energy_diffusion_C(pdt);
  clean_param_C(&three);
  T2rhoe_C();

  /* Prep for sources step - all workers must have a common average gravity and exospheric temperature */
  get_gavg_Tninf_C(&gavg,&Tninf);

  /* Sources substep and finalize solution for this time step */
  source_loss_allparams_C(pdt,pt,pymd,pUTsec,&f107a,&f107,pfirst,&gavg,&Tninf);
  clean_param_C(&three); clean_param_C(&two); clean_param_C(&one);

  // Fix electron veloc???
  return;
}
