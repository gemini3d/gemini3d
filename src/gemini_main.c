// MAIN PROGRAM FOR GEMINI3D

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

#include "gemini3d.h"

int main(int argc, char **argv) {

struct params s;

int ierr = MPI_Init(&argc, &argv);

if (argc < 2) {
  MPI_Finalize();
  fprintf(stderr, "Gemini3D: please give simulation output directory e.g. ~/data/my_sim\n");
  return 1;
}

int L = strlen(argv[1]);
if(L > LMAX) {
  MPI_Finalize();
  fprintf(stderr, "Gemini3D simulation output directory: path length > %d\n", LMAX);
  return 1;
}
L++; // for null terminator

if(strcmp(argv[1], "-") == 0) {
  MPI_Finalize();
  fprintf(stderr, "Gemini3d: first argument must be output directory\n");
  return 1;
}

strcpy(s.out_dir, argv[1]);

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
    return 0;
  }
  if (strcmp(argv[i], "-manual_grid") == 0) {
    if (argc < i+1) {
      MPI_Finalize();
      fprintf(stderr, "-manual_grid lid2in lid3in\n");
      return 1;
    }
    lid2in = atoi(argv[i]);
    lid3in = atoi(argv[i+1]);
  }
}

gemini_main(&s, &lid2in, &lid3in);

ierr = MPI_Finalize();

if (ierr != 0) return 1;

return 0;
}
