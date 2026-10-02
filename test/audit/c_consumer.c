/* A C17 client, compiled as C rather than C++. Apache-2.0. */
#include "gemini3d.h"
#include <stddef.h>
#include <stdio.h>
extern size_t audit_params_size(void);
int main(void) {
  if (sizeof(struct params) != audit_params_size()) return 1;
  void *cfg = NULL;
  gemini_cfg_alloc_C(&cfg);
  if (!cfg) return 2;
  struct guarded_int { int before, value, after; } bg = {12345, -1, 54321}, pert = bg;
  double dtbg = -1, dtpert = -1;
  get_config_vars_C(&cfg, &bg.value, &pert.value, &dtbg, &dtpert);
  if (bg.before != 12345 || bg.after != 54321 || pert.before != 12345 || pert.after != 54321) return 3;
  if (bg.value != 0 || pert.value != 0 || dtbg != 900 || dtpert != 0) return 4;
  int ymd[3] = {2024, 2, 28};
  double ut = 86399.5, step = 1.0;
  dateinc_C(&step, ymd, &ut);
  if (ymd[0] != 2024 || ymd[1] != 2 || ymd[2] != 29 || ut != 0.5) return 5;
  int species = 0;
  get_species_size_C(&species);
  if (species != 7) return 6;
  gemini_cfg_dealloc_C(&cfg);
  if (cfg) return 7;
  gemini_cfg_dealloc_C(&cfg);
  puts("C17 client layout, integer flags, calendar and lifetime passed");
  return 0;
}
