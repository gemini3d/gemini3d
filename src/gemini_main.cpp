// MAIN PROGRAM FOR GEMINI3D

#include <iostream>
#include <filesystem>
#include <vector>
#include <sstream>

#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "gemini3d.h"
#include "iniparser.h"
#include "pathlib.hpp"

namespace fs = std::filesystem;

int main(int argc, char **argv) {

struct params s;

int ierr = MPI_Init(&argc, &argv);

// CLI
if (argc < 2) {
  std::cerr << "Gemini3D: please give simulation output directory e.g. ~/data/my_sim" << std::endl;
  return EXIT_FAILURE;
}

// simulation directory
std::string out_dir = expanduser(argv[1]);
if(out_dir.size() > LMAX) {
  std::cerr << "Gemini3D simulation output directory: path length > " << LMAX << std::endl;
  return EXIT_FAILURE;
}


if(! fs::is_directory(out_dir)) {
  std::cerr << "Gemini3D simulation output directory does not exist: " << out_dir << std::endl;
  return EXIT_FAILURE;
}

// Read gemini_config.ini
std::string ini_file = out_dir;
ini_file.append("/inputs/gemini_config.ini");

dictionary  *ini;
int b,i ;
double d;
const char *txt;

ini = iniparser_load(ini_file.c_str());
if (ini==NULL) {
    std::cerr << "gemini3d_ini: cannot parse file: " << ini_file << std::endl;
    return EXIT_FAILURE;
}

std::string ini_str, t_str;
std::vector<int> ymd;

ini_str = iniparser_getstring(ini, "base:ymd", "");
if(ini_str.empty()) {
  std::cerr << "gemini3d_ini: base:ymd not found in " << ini_file << std::endl;
  return EXIT_FAILURE;
}
std::stringstream sini(ini_str);

while(std::getline(sini, t_str, ',')) ymd.push_back(stoi(t_str));
if(ymd.size() != 3) {
  std::cerr << "gemini3d_ini: base:ymd must have 3 elements: " << ini_str << std::endl;
  return EXIT_FAILURE;
}


iniparser_freedict(ini);  // close the file

// Prepare Gemini3D struct
std::strncpy(s.out_dir, out_dir.c_str(), LMAX);
s.fortran_cli = false;
s.debug = false;
s.dryrun = false;
int lid2in = -1, lid3in = -1;

for (int i = 2; i < argc; i++) {
  if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "-debug") == 0) s.debug = true;
  if (strcmp(argv[i], "-dryrun") == 0) s.dryrun = true;
  if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
    MPI_Finalize();
    help_gemini_bin();
    return EXIT_SUCCESS;
  }
  if (strcmp(argv[i], "-manual_grid") == 0) {
    if (argc < i+1) {
      MPI_Finalize();
      std::cerr << "-manual_grid lid2in lid3in" << std::endl;
      return EXIT_FAILURE;
    }
    lid2in = atoi(argv[i]);
    lid3in = atoi(argv[i+1]);
  }
}

gemini_main(&s, &lid2in, &lid3in);

ierr = MPI_Finalize();

if (ierr != 0) return EXIT_FAILURE;

return EXIT_SUCCESS;
}


// top-level module calls for gemini simulation
void gemini_main(struct params* ps, int* plid2in, int* plid3in){
  int ierr;
  int lx1,lx2,lx3;
  int lx2all,lx3all;
  int lsp;
  double UTsec;
  int[3] ymd; 
  double** fluidvars, fluidauxvars, electrovars;    // pointers modifiable by fortran
  double t=0, dt=1e-6;
  double tout, tneuBG, tglowout, tdur, tmilestone;
  double tstart, tfin;
  int it, iupdate;
  int flagoutput; 
  double v2grid,v3grid;
  bool first,flagneuBG
  int flagdneu;
  double dtneu,dtneuBG

  /* Basic setup */
  mpisetup();   // organize mpi workers
  cli_config_gridsize(ps,plid2in,plid3in);    // handling of input data, create internal fortran type with parameters for run
  get_fullgrid_size_C(&lx1,&lx2all,%lx3all)
  init_procgrid(&lx2all,&lx3all,plid2in,plid3in);    // compute process grid for this run
  get_config_vars_C(&flagneuBG,&flagdneu,&dtneuBG,&dtneu)

  /* Get input grid from file */
  read_grid_C();    // read the input grid from file, storage as fortran module object

  /* Main needs to know the grid sizes and species numbers */
  get_subgrid_size(&lx1,&lx2,&lx3);   // once grid is input we need to know the sizes
  get_species_size(&lsp);

  /* Allocate memory and get pointers to blocks of data */
  gemini_alloc(fluidvars,fluidauxvars,electrovars);
  outdir_fullgridvaralloc(&lx1,&lx2all,&lx3all);

  /* initialize state variables from input file */ 
  get_initial_state(&UTsec,&ymd[0],&tdur);
  set_start_values(&it,&t,&tout,&tglowout,&tneuBG);

  /* initialize other file input data */
  fprintf(" Initializing electric field input data...\n");
  init_Efieldinput(&dt,&t,&ymd[0],&UTsec);
  pot2perpfield_C();
  BGfield_Lagrangian(&v2grid,&v3grid);
  fprintf(" Initialize precipitation input data...\n");
  init_precipinput(&dt,&t,&ymd[0],&UTsec);
  fprintf(" Initialize neutral background and input files...\n");
  msisinit_C();
  init_neutralBG_C(&dt,&t,&ymd[0],&UTsec,&v2grid,&v3grid);
  init_neutralperturb_C(&dt,&ymd[0],&UTsec);

  /* Compute initial drift velocity */
  get_initial_drifts();

  /* Control console printing for, actually superfluous FIXME */
  set_update_cadence(&iupdate);

  while(t<tdur){
    dt_select_C();

    // neutral data
    if (it~=1 .and. flagneuBG .and. t>tneuBG){
      neutral_atmos_winds_C(&ymd[0],&UTsec);
      neutral_atmos_wind_update_C(&v2grid,&v3grid);
      tneuBG+=dtneuBG;
      fprintf(" Computed neutral background...\n");
    }
    if (flagdneu==1){
      neutral_perturb_C(&dt,&t,&ymd[0],&UTsec,&v2grid,&v3grid);
      fprintf(" Computed neutral perturbations...\n");
    }

    // call electrodynamics solution
    electrodynamics_C(&it,&t,&dt,&ymd[0],&UTsec);
    fprintf(" Computed electrodynamics solutions...\n");

    // advance the fluid state variables
    first=it==1;
    fluid_adv(&t,&dt,&ymd[0],&UTsec,&first,&lsp);
    fprint(" Computed fluid update...\n");

    check_finite_output(&t);
    it=+1; t+=dt; 
    fprintf(" Time step finished:  %d %d %d %f",ymd[0],ymd[1],ymd[2],UTsec);
    check_dryrun();
    check_fileoutput(&t,&tout,&tglowout,&tmilestone,&flagoutput,&ymd[0],&UTsec); 
  }

  /* Call deallocation procedures */
  gemini_dealloc(fluidvars,fluidauxvars,electrovars);
  clear_neutralBG_C();
  clear_dneu_C();

  return;
}


void fluid_adv(double* pt, double* pdt, int* ymd, double* pUTsec, bool* pfirst){
  double f107,f107a;
  double gavg,Tninf;
  int one=1,two=2,three=3;    // silly but I need some way to pass these ints by reference to fortran...

  /* Set up variables for the time step */ 
  get_solar_indices_C(&f107,&f107a);   // FIXME: do we really need to return the indices???
  v12rhov1();
  T2rhoe();

  /* Advection substep */
  halo_interface_vels_allspec_C(plsp);
  interface_vels_allspec(plsp);
  set_global_boundaries_allspec_C(plsp);
  halo_allparams_C();
  sweep3_allparams_C(pdt);
  sweep1_allparams_C(pdt);
  halo_allparams_C();
  sweep2_allparams_C(pdt);
  rho12v1_C();
  clean_param_C(&one);
  claen_param_C(&two);

  /* Compression substep */
  VNRicht_artvisc_C();
  RK2_prep_mpi_allspec_C();
  compression(pdt);
  rhoe2T();
  clean_param(&three);

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
