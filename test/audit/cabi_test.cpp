#include "gemini3d.h"
#include <cstddef>
#include <iostream>
#include <type_traits>
extern "C" std::size_t audit_params_size();
int main(int argc,char**) {
  static_assert(std::is_standard_layout<params>::value);
  if(audit_params_size()!=sizeof(params)) return 1;
  void *cfg=nullptr;
  if(argc>1) {
    params p{};
    gemini_cfg_alloc_C(&cfg);
    read_config_in_C(&p,&cfg); // Unsupported path must terminate with an explicit error.
    return 0;
  }
  for(int i=0;i<100;++i) {
    gemini_cfg_alloc_C(&cfg);
    struct Guard {int left;int value;int right;};
    Guard a{0x12ab34cd,-1,0x34cd12ab},b=a;
    double dtbg=-1,dt=-1;
    get_config_vars_C(&cfg,&a.value,&b.value,&dtbg,&dt);
    if(a.left!=0x12ab34cd||a.right!=0x34cd12ab||b.left!=a.left||b.right!=a.right) return 2;
    if(a.value!=0||b.value!=0||dtbg!=900||dt!=0) return 3;
    gemini_cfg_dealloc_C(&cfg);
    if(cfg!=nullptr)return 4;
    gemini_cfg_dealloc_C(&cfg); // Idempotent null deallocation.
  }
  std::cout << "C/Fortran ABI and handle checks passed\n";
}
