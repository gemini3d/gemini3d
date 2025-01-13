// MAIN PROGRAM FOR GEMINI3D

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <exception>

#include <filesystem>

namespace fs = std::filesystem;

#include <mpi.h>

#include "gemini3d.h"
#include "ffilesystem.h"


int main(int argc, char **argv) {

  struct params s;
  int myid;
  int ierr = MPI_Init(&argc, &argv);
  if(ierr){
    std::cerr << "MPI_Init failed\n";
    return EXIT_FAILURE;
  }

  // CLI
  if (argc < 2) {
    help_gemini_bin();
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::string_view a1(argv[1]);

  if (a1 == "-h" || a1 == "-help") {
    help_gemini_bin();
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  // simulation directory
  std::string out_dir(fs_expanduser(a1));

  if( !fs_is_dir(out_dir)){
    std::cerr << "Gemini3D simulation output directory does not exist: " << out_dir << "\n";
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // we don't have a C++ parser for Fortran namelist files,
  // so read the namelist file directly in Fortran as usual.
  s.fortran_nml = true;

  // Prepare Gemini3D struct
  std::strcpy(s.out_dir, out_dir.data());

  s.fortran_cli = false;
  s.debug = false;
  s.dryrun = false;
  int lid2in = -1, lid3in = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  for (int i = 2; i < argc; i++) {
    std::string_view arg(argv[i]);
    if (arg == "-d" || arg == "-debug")
      s.debug = true;

    if (arg == "-dryrun")
      s.dryrun = true;

    if (arg == "-h" || arg == "-help") {
      help_gemini_bin();
      MPI_Finalize();
      return EXIT_SUCCESS;
    }
    if (arg == "-manual_grid") {
      if (argc < i+1) {
        MPI_Finalize();
        std::cerr << "-manual_grid lid2in lid3in\n";
        return EXIT_FAILURE;
      }
      lid2in = atoi(argv[i]);
      lid3in = atoi(argv[i+1]);
    }
  }

  gemini_main(&s, &lid2in, &lid3in);

  if(MPI_Finalize()){
    std::cerr << "MPI_Finalize failed\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
