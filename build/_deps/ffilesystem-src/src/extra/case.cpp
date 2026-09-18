#include "ffilesystem.h"

#include <string_view>
#include <string>


namespace {
struct fs_remove_on_exit {
  std::string path;

  ~fs_remove_on_exit()
  {
    if(!path.empty())
      fs_remove(path);
  }
};

}

bool fs_is_case_sensitive(std::string_view path)
{
  // check if the filesystem is case sensitive by creating a temporary file
  // with a mixed-case name and checking if the same name with different case
  // resolves to the same file.

  if(path.empty() || !fs_is_dir(path))
    return false;

  constexpr std::string::size_type N = 16; // arbitrary

  const std::string rname = fs_generate_random_alphanumeric_string(N);
  if(rname.empty())
    return false;

  const std::string mixed_name = ".FfilesystemCaseSensitiveCheck_" + rname;
  const std::string mixed_path = fs_join(path, mixed_name);

  if(!fs_touch(mixed_path))
    return false;

  const fs_remove_on_exit cleanup{mixed_path};

  std::string lower_name = mixed_name;
  fs_ascii_lower(lower_name);

  const std::string lower_path = fs_join(path, lower_name);

  // Case-insensitive filesystems resolve both names to the same file.
  return fs_equivalent(mixed_path, lower_path);
}
