// MAIN PROGRAM FOR GEMINI3D

#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#else
#error "No C++ filesystem support"
#endif

#include <mpi.h>

#include "gemini3d.h"
#include "ffilesystem.h"

int main(int argc, char **argv) {

  struct params s;
  int myid;
  int ierr = MPI_Init(&argc, &argv);

  // CLI
  if (argc < 2) {
    help_gemini_bin();
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
    help_gemini_bin();
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  // simulation directory
  char odir[4096];
  fs_expanduser(argv[1], odir, 4096);
  fs::path out_dir(odir);

  if(! fs::is_directory(out_dir)) {
    std::cerr << "Gemini3D simulation output directory does not exist: " << out_dir << std::endl;
    return EXIT_FAILURE;
  }

  // Read gemini_config.ini, if it exists
  auto ini_file = out_dir / "inputs/gemini_config.ini";
  if(fs::is_regular_file(ini_file)) {

    dictionary  *ini;

    // TODO: use libsc ini parser
    ini = iniparser_load(ini_file.string().c_str());
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

    s.fortran_nml = false;
  }
  else {
    s.fortran_nml = true;
  }

  // Prepare Gemini3D struct
  std::strcpy(s.out_dir, out_dir.generic_string().c_str());

  s.fortran_cli = false;
  s.debug = false;
  s.dryrun = false;
  int lid2in = -1, lid3in = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  for (int i = 2; i < argc; i++) {
    if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "-debug") == 0) s.debug = true;
    if (strcmp(argv[i], "-dryrun") == 0) s.dryrun = true;
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
      help_gemini_bin();
      MPI_Finalize();
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
