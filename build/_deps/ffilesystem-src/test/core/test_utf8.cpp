#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"intl"_test = [] {
const std::string smiley = "😀";
const std::string wink = "😉";
const std::string hello = "你好";

expect(!fs_canonical(".", true).empty());
expect(!fs_canonical("./", true).empty());

expect(eq(fs_file_name("./" + smiley), smiley));
expect(eq(fs_file_name("./" + wink), wink));
expect(eq(fs_file_name("./" + hello), hello));

// This is shaky on macOS. Better to leave it off for now.

// EXPECT_EQ(fs_canonical(smiley, false, false), smiley);
// EXPECT_EQ(fs_canonical(wink, false, false), wink);
// EXPECT_EQ(fs_canonical(hello, false, false), hello);


};
}
