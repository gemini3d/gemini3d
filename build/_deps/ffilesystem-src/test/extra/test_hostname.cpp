#include "ffilesystem.h"

#include <vector>
#include <string>

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"Hostname"_test = [] {
  std::string s = fs_hostname();

  expect(!s.empty());

  expect(s.length() != fs_get_max_path()) << "hostname length is equal to max path length";
};

"MaxComponent"_test = [] {
  expect(fs_max_component("/") >= 1);
};


"Shell"_test = [] {
  std::string s = fs_get_shell();
  if (!s.empty())
    expect(s.length() != fs_get_max_path()) << "shell has blank space";
};


"Terminal"_test = [] {
  std::string s = fs_get_terminal();
  if (!s.empty())
    expect(s.length() != fs_get_max_path()) << "terminal has blank space";
};
}
