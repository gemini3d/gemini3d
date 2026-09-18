#include "ffilesystem.h"

#include <array>
#include <string_view>

#include <boost/ut.hpp>

int main() {
	using namespace boost::ut;

	"drop_slash"_test = [] {
		struct case_t {
			std::string_view input;
			std::string_view expected;
		};

		constexpr std::array<case_t, 9> cases{{
			{"", ""},
			{"a", "a"},
			{"a/", "a"},
			{"a/b", "a/b"},
			{"a/b/", "a/b"},
			{"////", "/"},
			{"a////b", "a/b"},
			{"a//b//", "a/b"},
			{"/", "/"},
		}};

		for (const auto& test_case : cases) {
			expect(eq(fs_drop_slash(test_case.input), test_case.expected));
		}

		if (fs_is_windows()) {
			constexpr std::array<case_t, 3> windows_cases{{
				{"c:/", "c:/"},
				{"c:///", "c:/"},
				{"c:/a/b//", "c:/a/b"},
			}};

			for (const auto& test_case : windows_cases) {
				expect(eq(fs_drop_slash(test_case.input), test_case.expected));
			}

			if (fs_win32_long_paths_enabled()) {
				expect(eq(fs_drop_slash(R"(\\?\C:/)"), std::string_view{R"(\\?\C:/)"}));
			}
		}
	};
}
