#include <boost/ut.hpp>

#include "ffilesystem.h"

#include <cctype> // for std::isalnum

int main() {
  using namespace boost::ut;

  "alphanumeric_string"_test = [] {
  const std::string s = fs_generate_random_alphanumeric_string(16);
  expect(eq(s.size(), std::size_t{16}));
  for (char c : s) {
    expect(std::isalnum(static_cast<unsigned char>(c)) != 0);
  }
  };
}
