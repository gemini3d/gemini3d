#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"root"_test = [] {

expect(eq(fs_root(""), std::string{""}));
expect(eq(fs_root("a/b"), std::string{""}));
expect(eq(fs_root("./a/b"), std::string{""}));
expect(eq(fs_root("../a/b"), std::string{""}));

if(fs_is_windows()){
expect(eq(fs_root("c:"), std::string{"c:"}));
expect(eq(fs_root("c:/a/b"), std::string{"c:/"}));
expect(eq(fs_root("/etc"), std::string{"/"}));
expect(eq(fs_root("\\etc"), std::string{"\\"}));
expect(eq(fs_root("c:\\"), std::string{"c:\\"}));
expect(eq(fs_root("c:/"), std::string{"c:/"}));
expect(eq(fs_root("\\"), std::string{"\\"}));
} else {
expect(eq(fs_root("/a/b"), std::string{"/"}));
expect(eq(fs_root("c:/etc"), std::string{""}));
}

};

"root_name"_test = [] {
expect(eq(fs_root_name(""), std::string{""}));
expect(eq(fs_root_name("a/b"), std::string{""}));
expect(eq(fs_root_name("./a/b"), std::string{""}));
expect(eq(fs_root_name("../a/b"), std::string{""}));

if (fs_is_windows()){
expect(eq(fs_root_name("c:/a/b"), std::string{"c:"}));
expect(eq(fs_root_name("/etc"), std::string{""}));
expect(eq(fs_root_name(R"(C:\)"), std::string{"C:"}));

const std::string drive_prefixed = "C:/must-not-be-read";
const std::string_view truncated_drive(drive_prefixed.data(), 1);
expect(eq(truncated_drive, std::string_view{"C"}) >> fatal);
expect(eq(fs_root_name(truncated_drive), std::string{""}))
  << "fs_root_name() read past string_view length for a non-null-terminated path";
} else {
expect(eq(fs_root_name("/a/b"), std::string{""}));
expect(eq(fs_root_name("c:/etc"), std::string{""}));
}

};
}
