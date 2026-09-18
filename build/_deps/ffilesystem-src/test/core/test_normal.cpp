#include "ffilesystem.h"
#include <array>
#include <string>

#include <boost/ut.hpp>

int main() {
    using namespace boost::ut;

    struct case_t {
        std::string inp;
        std::string expected;
    };

    const std::array<case_t, 48> cases{{
            {"", "."},
            {"/", "/"},
            {"//", "/"},
            {"/////", "/"},
            {".", "."},
            {"./", "."},
            {"./.", "."},
            {"..", ".."},
            {"../", ".."},
            {"a/..", "."},
            {"../..", "../.."},
            {"a/b/..", "a"},
            {"a/b/../..", "."},
            {"a/b/../../..", ".."},
            {"/a", "/a"},
            {"/a/", "/a"},
            {"/a/.", "/a"},
            {"/a/..", "/"},
            {"/a/b/..", "/a"},
            {"a", "a"},
            {".a", ".a"},
            {"a.", "a."},
            {"a./", "a."},
            {"a/b", "a/b"},
            {"..a", "..a"},
            {"a..", "a.."},
            {"a../", "a.."},
            {"a/", "a"},
            {"a//", "a"},
            {"./a", "a"},
            {"./a/", "a"},
            {"./a/.", "a"},
            {"../a", "../a"},
            {"../a/b/..", "../a"},
            {"a/b/", "a/b"},
            {"a/b/.", "a/b"},
            {"a/b/..", "a"},
            {"a/b/../", "a"},
            {"a/b/../c", "a/c"},
            {"a/b/../c/d", "a/c/d"},
            {"a/b/../../c/d", "c/d"},
            {"././a/./b/././c/./.", "a/b/c"},
            {"a/b/../../c/../..", ".."},
            {"a/b/../../../c/../..", "../.."},
            {"a/./b/..", "a"},
            {"a/.///b/../", "a"},
            {"/a/../..", "/"},
            {"/a/../../../../", "/"},
    }};

    "normalize"_test = [cases] {
        for (const auto& test_case : cases) {
            expect(eq(fs_normal(test_case.inp), test_case.expected));
        }
    };

#if defined(_WIN32)
    const std::array<case_t, 7> windows_cases{{
            {R"(\/\///\/)", "/"},
            {R"(a/b/..\//..///\/../c\\/)", "../c"},
            {R"(..a/b/..\//..///\/../c\\/)", "../c"},
            {R"(..\)", ".."},
            {R"(c:\)", "c:/"},
            {R"(c:\\)", "c:/"},
            {R"(c:\a/b/../)", "c:/a"},
    }};

    "normalize_windows"_test = [windows_cases] {
        for (const auto& test_case : windows_cases) {
            expect(eq(fs_normal(test_case.inp), test_case.expected));
        }
    };
#endif
// some tests from https://github.com/gulrak/filesystem/blob/b1982f06c84f08a99fb90bac43c2d03712efe921/test/filesystem_test.cpp#L950
}
