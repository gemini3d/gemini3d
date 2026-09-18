#include "ffilesystem.h"
#include <iostream>

#include <boost/ut.hpp>


int main() {
using namespace boost::ut;

"proximate"_test = [] {

expect(eq(fs_proximate_to("", ""), std::string{""}));
expect(eq(fs_proximate_to("Hello", "Hello"), std::string{"."}));
expect(eq(fs_proximate_to("Hello", "Hello/"), std::string{"."}));
expect(eq(fs_proximate_to("a/b", "a/b"), std::string{"."}));
expect(eq(fs_proximate_to("a///b", "a/b"), std::string{"."}));
expect(eq(fs_proximate_to("a/b", "a///b"), std::string{"."}));
expect(eq(fs_proximate_to("a/b", "a/b/"), std::string{"."}));
expect(eq(fs_proximate_to("a/b/", "a/b"), std::string{"."}));
expect(eq(fs_proximate_to("a/b", "a"), std::string{".."}));
expect(eq(fs_proximate_to("a/b", "a/"), std::string{"../"}));
expect(eq(fs_proximate_to("a", "a/b/"), std::string{"b/"}));
expect(eq(fs_proximate_to("a", "a/b/."), std::string{"b/"}));
expect(eq(fs_proximate_to("a", "a/b/.."), std::string{"."}));
expect(eq(fs_proximate_to("a/b/c/d", "a/b"), std::string{"../.."}));
expect(eq(fs_proximate_to("a/b/c/d", "a/b/"), std::string{"../../"}));
expect(eq(fs_proximate_to("./a/b", "./a/c"), std::string{"../c"}));

if (fs_is_windows()){
auto e = fs_getenv("SYSTEMDRIVE");
expect(e.has_value() >> fatal) << "Failed to get SYSTEMDRIVE environment variable";
std::string c = e.value();

expect(eq(fs_proximate_to(c+"/", c+"/a/b"), std::string{"a/b"}));
expect(eq(fs_proximate_to(c+"/a/b", c+"/a/b"), std::string{"."}));
expect(eq(fs_proximate_to(c+"/a/b", c+"/a"), std::string{".."}));
// {c+"/a", "b", "b"}  //  ambiguous with Clang/Flang ARM MinGW <filesystem>

// NOTE: on Windows, if a path is real, finalPath is used, which makes drive letters upper case.
// we didn't use a random drive letter because a removable drive at that letter can make the test fail spuriously.

} else {
expect(eq(fs_proximate_to("", "a"), std::string{"a"}));
expect(eq(fs_proximate_to("/", "/"), std::string{"."}));
expect(eq(fs_proximate_to("Hello", "Hello"), std::string{"."}));
expect(eq(fs_proximate_to("Hello", "Hello/"), std::string{"."}));
expect(eq(fs_proximate_to("/dev/null", "/dev/null"), std::string{"."}));
expect(eq(fs_proximate_to("a/b", "c"), std::string{"../../c"}));
expect(eq(fs_proximate_to("c", "a/b"), std::string{"../a/b"}));
expect(eq(fs_proximate_to("a/b", "a/c"), std::string{"../c"}));
}
// NOTE: use relative non-existing paths, as on macOS AppleClang, the <filesystem> gives incorrect results on non-existing absolute paths,
// Which don't make sense anyway.

};
}
