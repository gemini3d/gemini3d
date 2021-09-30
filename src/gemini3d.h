#ifndef GEMINI3D_H
#define GEMINI3D_H

enum { LMAX = 1000 };

struct params {
  // order and lengths must match in Fortran and C
  bool fortran_cli;
  bool debug;
  bool dryrun;
  char out_dir[LMAX];
};


extern void gemini_main(struct params, int*, int*);
#endif
