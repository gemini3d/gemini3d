#include <ISO_Fortran_binding.h>

int realpath_cfi(const CFI_cdesc_t *in_desc, CFI_cdesc_t *out_desc) {
    if(!in_desc || !out_desc || !in_desc->base_addr || !out_desc->base_addr)
      return -1;

    if(in_desc->rank != 0 || out_desc->rank != 0)
      return -2;

    if(in_desc->type != CFI_type_char || out_desc->type != CFI_type_char)
      return -3;

    return 0;
}
