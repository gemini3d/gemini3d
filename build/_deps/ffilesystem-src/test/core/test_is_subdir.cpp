#include "ffilesystem.h"
#include <array>
#include <string>

#include <boost/ut.hpp>

int main() {
    using namespace boost::ut;

    struct case_t {
        std::string base;
        std::string sub;
        bool expected;
    };

    const std::array<case_t, 11> cases{{
            {"a/b/c", "a/b", true},
            {"a/b/c", "a/b/", true},
            {"a/b/c", "a", true},
            {"a/b/", "a/b", false},
            {"a/b", "a/b", false},
            {"/a/b", "a/b", false},
            {"a/b", "/a/b", false},
            {"a/b", "a/b/", false},
            {"a/b", "c", false},
            {"b", "a/b", false},
            {"c:/a", "c:/", true},
    }};

    "is_subdir"_test = [cases] {
        for (const auto& test_case : cases) {
            expect(eq(fs_is_subdir(test_case.base, test_case.sub), test_case.expected));
        }
    };
}
