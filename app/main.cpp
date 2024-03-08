// MAIN PROGRAM FOR GEMINI3D

#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>
#include <exception>

#include <filesystem>
static_assert(__cpp_lib_filesystem, "C++17 <filesystem> required");

namespace fs = std::filesystem;

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
  fs::path out_dir(Ffs::expanduser(argv[1]));

  if(! fs::is_directory(out_dir))
    throw std::runtime_error("Gemini3D simulation output directory does not exist: " + out_dir.string());

  // we don't have a C++ parser for Fortran namelist files,
  // so read the namelist file directly in Fortran as usual.
  s.fortran_nml = true;

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
        throw std::runtime_error("-manual_grid lid2in lid3in");
      }
      lid2in = atoi(argv[i]);
      lid3in = atoi(argv[i+1]);
    }
  }

  gemini_main(&s, &lid2in, &lid3in);

  if(MPI_Finalize())
    throw std::runtime_error("MPI_Finalize failed");

  return EXIT_SUCCESS;
}
