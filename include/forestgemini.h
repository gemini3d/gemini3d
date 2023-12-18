#ifndef FORESTGEMINI_H
#define FORESTGEMINI_H

#ifdef __cplusplus
extern "C" {
#endif

extern void grid_from_extents_C(double*,double*,double*,double*,double*,int*,int*,int*,void**);
extern void gemini_grid_alloc_C(double*,double*,double*,int*,int*,int*,int*,void**);
extern void interp_file2subgrid_C(void**, int*, void**, double**, double**);
extern void clean_param_after_regrid_C(int*, int*, void**, double**,void**);
extern void checkE1_C(double**, double**, double**, int*);
extern void forceZOH_all_C(double**);
extern void permute_fluidvars_C(double**);
extern void ipermute_fluidvars_C(double**);
extern void tag4refine_C(int*,void**,double**,double**,double**,void**,int*,bool*);
extern void tag4coarsening_C(int*,void**,double**,double**,double**,void**,bool*);
extern void get_grid_magcoords_C(int*,void**,double**,double**,double**);
extern void get_grid_magcoordsi_C(int*,void**,double**,double**,double**);
extern void get_locationsi_C(void**,bool*,double*,double*,double*,double**,double**,double**,double**,int*,int*);
extern void get_datainow_ptr_C(void**,double**);
extern void set_datainow_C(void**);
extern void swap_statevars_C(double**, double**);

#ifdef __cplusplus
}
#endif

#endif
