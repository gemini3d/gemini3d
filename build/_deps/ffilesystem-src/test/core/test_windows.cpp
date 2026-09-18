#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

  "ToWide"_test = [] {
    std::string_view s = "hello";
    std::wstring w = fs_win32_to_wide(s);

    expect(w.size() == s.size());
    expect(w == L"hello");
  };

  "ToNarrow"_test = [] {
    std::wstring_view w = L"hello";
    std::string n = fs_win32_to_narrow(w);

    expect(n.size() == w.size());
    expect(n == "hello");
  };
}
