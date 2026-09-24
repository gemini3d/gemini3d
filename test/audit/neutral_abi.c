#include "gemini3d.h"
#include <stdlib.h>

extern void audit_neutral_fixture(void **, void **, void **, const int *);
extern void audit_neutral_check(void **, const double *, const double *);
extern void audit_neutral_release(void **);

int main(int argc, char **argv) {
  int version = argc > 1 ? atoi(argv[1]) : 0;
  void *cfg = NULL, *grid = NULL, *work = NULL;
  int xtype = 1, ymd[3] = {2013, 2, 20};
  double ut = 18000, v2 = 12, v3 = -7;
  audit_neutral_fixture(&cfg, &grid, &work, &version);
  setv2v3_C(&v2, &v3);
  msisinit_C(&cfg);
  neutral_atmos_winds_C(&cfg, &xtype, &grid, ymd, &ut, &work);
  audit_neutral_check(&work, &v2, &v3);
  v2 = -18;
  v3 = 23;
  setv2v3_C(&v2, &v3);
  for (int i = 0; i < 2; ++i) {
    neutral_atmos_wind_update_C(&work);
    audit_neutral_check(&work, &v2, &v3);
  }
  audit_neutral_release(&work);
  return 0;
}
