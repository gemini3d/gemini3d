// MAIN PROGRAM FOR GEMINI3D

#include <iostream>
#include <filesystem>
#include <vector>
#include <sstream>

#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include <stdio.h>

#include "gemini3d.h"
#include "iniparser.h"
#include "filesystem.h"

namespace fs = std::filesystem;

int main(int argc, char **argv) {

  struct params s;
  int myid;
  int ierr = MPI_Init(&argc, &argv);

  // CLI
  if (argc < 2) {
    std::cerr << "Gemini3D: please give simulation output directory e.g. ~/data/my_sim" << std::endl;
    return EXIT_FAILURE;
  }

  // simulation directory
  char odir[4096];
  expanduser(argv[1], odir);
  std::string out_dir(odir);
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

  //std::strncpy(s.out_dir,"~/simulations/raid/CI/GDI_periodic_lowres_C_optim/",LMAX);
  s.fortran_cli = false;
  s.debug = false;
  s.dryrun = false;
  int lid2in = -1, lid3in = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

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
