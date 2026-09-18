#include "ffilesystem.h"
#include <array>
#include <string>

#include <boost/ut.hpp>

int main() {
    using namespace boost::ut;

    struct case_t {
        std::string base;
        std::string pre;
        bool expected;
    };

    const std::array<case_t, 12> cases{{
            {"a/b//c", "a/b", false},
            {"a/b/c", "a/b/", false},
            {"a/b/c", "a", false},
            {"a/b", "a/b", true},
            {"a/b/", "a/b", true},
            {"/a/b", "a/b", false},
            {"a/b", "/a/b", false},
            {"a/b", "a/b/", true},
            {"a/b", "c", false},
            {"b", "a/b", false},
            {"c:/a", "c:/", false},
            {"c:/", "c:/a", true},
    }};

    "is_prefix"_test = [cases] {
        for (const auto& test_case : cases) {
            expect(eq(fs_is_prefix(test_case.base, test_case.pre), test_case.expected));
        }
    };
}
