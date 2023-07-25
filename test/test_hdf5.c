#include "hdf5.h"
#include "hdf5_hl.h"

#include <stdlib.h>
#include <stdio.h>


int main(void){

/*
check that repeated calls to h5open_f() do not cause problems as per docs
not calling h5open_f() at all makes failures as library isn't initialized
unlike C HDF5, Fortran HDF5 does not auto-initialize.
*/

hid_t fid;
int rank = 1;
hsize_t dims[1] = {1};

float buf[1] = {42.0};

unsigned int major, minor, release;

if(H5get_libversion(&major, &minor, &release) < 0){
  fprintf(stderr, "ERROR:hdf5_standalone_C:H5get_libversion: could not get HDF5 library version\n");
  return EXIT_FAILURE;
}
printf("hdf5_standalone_C: HDF5 library version %d.%d.%d\n", major, minor, release);


if ( (fid = H5Fcreate("test_hdf5_standalone_C.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) == H5I_INVALID_HID){
  fprintf(stderr, "ERROR:hdf5_standalone_C: could not create file\n");
  return EXIT_FAILURE;
}
printf("hdf5_standalone_C: created file\n");

if (H5LTmake_dataset_float(fid, "A", rank, dims, buf) < 0){
  fprintf(stderr, "ERROR:hdf5_standalone_C: could not create dataset\n");
  return EXIT_FAILURE;
}
printf("hdf5_standalone_C: created dataset\n");

if (H5Fclose(fid) < 0){
  fprintf(stderr, "ERROR:hdf5_standalone_C: could not close file\n");
  return EXIT_FAILURE;
}
printf("hdf5_standalone_C: closed file\n");

printf("OK: hdf5_standalone_C\n");

return EXIT_SUCCESS;

}
