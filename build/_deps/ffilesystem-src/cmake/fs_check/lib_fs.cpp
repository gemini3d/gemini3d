#include <cstring>
#include <filesystem>

#include "myfs.h"


bool has_filename(const char* path){
  std::filesystem::path p(path);
  return p.has_filename();
}
