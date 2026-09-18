#include <string>

#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

suite TestLocale = [] {
  "LocaleName"_test = [] {
    std::string loc = fs_get_locale_name();
    // macOS and MSVC have empty locale. Don't fail
    if(!fs_is_macos() && !fs_is_msvc()){
      expect(!loc.empty());
    }
  };
};
}
