#include <string>
#include <array>

#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
  using namespace boost::ut;

  struct case_t {
    std::string inp;
    std::string exp;
  };

  const std::array<case_t, 16> cases{{
      {"", "."},
      {"/", "/"},
      {".", "."},
      {"./", "."},
      {"..", "."},
      {"../", "."},
      {"a", "."},
      {"a/", "."},
      {"a/.", "a"},
      {"a/..", "a"},
      {"a/b", "a"},
      {"a/b/", "a"},
      {"a/b/c", "a/b"},
      {"ab/.parent", "ab"},
      {"ab/.parent.txt", "ab"},
      {"a/b/../.parent.txt", "a/b/.."},
  }};

  "parent"_test = [cases] {
    for (const auto& test_case : cases) {
      expect(eq(fs_parent(test_case.inp), test_case.exp));
    }
  };

if (!fs_is_windows()) {
  skip / "parent_windows"_test = [] {};
} else {
  const std::array<case_t, 4> windows_cases{{
      {"c:\\a\\b/../.parent.txt", "c:\\a\\b/.."},
      {"x:/", "x:/"},
      {"x:\\", "x:/"},
      {R"(\\?\C:\a\b/../.parent.txt)", R"(\\?\C:\a\b/..)"},
  }};

  "parent_windows"_test = [windows_cases] {
    for (const auto& test_case : windows_cases) {
      if (fs_backend() == "<filesystem>" && fs_win32_is_ext_path(test_case.inp)) {
        return;
      }
      expect(eq(fs_parent(test_case.inp), test_case.exp));
    }
  };
}

}
