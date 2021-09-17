// MAIN PROGRAM FOR GEMINI3D

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>

#include "gemini3d.h"

int main(int argc, char **argv) {

int ierr = MPI_Init(&argc, &argv);

if (argc < 2) {
  perror("please give simulation output directory e.g. ~/data/my_sim");
  return 1;
}

int Lpath_max = 10000;
int L = strlen(argv[1]);
if(L > Lpath_max) {
  // arbitrary maximum path length--10000 is not referred to anywhere else in the program.
  fprintf(stderr, "Gemini3D simulation output directory: path length > %d", Lpath_max);
  return 1;
}
L++; // for null terminator

char *out_dir = (char *)malloc(L);
if (out_dir != NULL) {
  strcpy(out_dir, argv[1]);
}

int lid2in = -1, lid3in = -1;
bool fortran_cli = false;

gemini_main(out_dir, &L, &lid2in, &lid3in, &fortran_cli);

ierr = MPI_Finalize();

if (ierr != 0) return 1;

return 0;
}
