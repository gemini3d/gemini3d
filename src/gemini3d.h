#ifndef GEMINI3D_H
#define GEMINI3D_H

#ifdef __cplusplus
extern "C" {
#endif

enum { LMAX = 1000 };

struct params {
  // order and lengths must match in Fortran and C
  // see gemini_main.f90 "cparams"
  bool fortran_cli;
  bool debug;
  bool dryrun;
  char out_dir[LMAX];
  // .ini [base]
  int ymd[3];
  float UTsec0;
  float tdur;
  float dtout;
  float activ[3];
  float tcfl;
  float Teinf;
  // .ini
};


extern void gemini_main(struct params *, int*, int*);

extern void help_gemini_bin();

#ifdef __cplusplus
}
#endif

#endif
