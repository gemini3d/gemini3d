#include "ffilesystem.h"

#include <boost/ut.hpp>

int main() {
using namespace boost::ut;

"stem"_test = [] {

expect(eq(fs_stem(""), std::string{""}));
expect(eq(fs_stem("stem.a.b"), std::string{"stem.a"}));
expect(eq(fs_stem("stem.a"), std::string{"stem"}));
expect(eq(fs_stem("stem"), std::string{"stem"}));
expect(eq(fs_stem(".stem"), std::string{".stem"}));
expect(eq(fs_stem("stem."), std::string{"stem"}));
expect(eq(fs_stem("stem.a."), std::string{"stem.a"}));
expect(eq(fs_stem("stem/a/b"), std::string{"b"}));
expect(eq(fs_stem("./.stem"), std::string{".stem"}));
expect(eq(fs_stem("../.stem"), std::string{".stem"}));
expect(eq(fs_stem(".stem.txt"), std::string{".stem"}));
expect(eq(fs_stem("a/.."), std::string{".."}));
expect(eq(fs_stem("a/../"), std::string{""}));
expect(eq(fs_stem("a/."), std::string{"."}));

expect(eq(fs_stem("日本語.日本語"), std::string{"日本語"}));
expect(eq(fs_stem("some space.txt"), std::string{"some space"}));
expect(eq(fs_stem("some space"), std::string{"some space"}));

if(fs_is_windows()){

expect(eq(fs_stem(R"(C:\a\ball.text)"), std::string{"ball"}));

if(fs_win32_long_paths_enabled()) {

expect(eq(fs_stem(R"(\\?\)"), std::string{""}));
expect(eq(fs_stem(R"(\\.\)"), std::string{""}));

expect(eq(fs_stem(R"(\\?\C:\)"), std::string{""}));
expect(eq(fs_stem(R"(\\.\C:\)"), std::string{""}));
expect(eq(fs_stem(R"(\\?\UNC\server\share)"), std::string{"share"}));
expect(eq(fs_stem(R"(\\?\UNC\server\share\日本語.txt)"), std::string{"日本語"}));

expect(eq(fs_stem(R"(\\server\share\some space here.txt)"), std::string{"some space here"}));
expect(eq(fs_stem(R"(\\?\C:\some space here.txt)"), std::string{"some space here"}));
}
}

};

"suffix"_test = [] {
expect(eq(fs_suffix(""), std::string{""}));
expect(eq(fs_suffix("a"), std::string{""}));
expect(eq(fs_suffix("a."), std::string{"."}));
expect(eq(fs_suffix("a.b"), std::string{".b"}));
expect(eq(fs_suffix("a.b.c"), std::string{".c"}));
expect(eq(fs_suffix("a/b.c"), std::string{".c"}));
expect(eq(fs_suffix("a/b.c.d"), std::string{".d"}));
expect(eq(fs_suffix("a/b/c.d"), std::string{".d"}));
expect(eq(fs_suffix("a/b/c.d.e"), std::string{".e"}));
expect(eq(fs_suffix("a/b/c.d/e"), std::string{""}));
expect(eq(fs_suffix(".a"), std::string{""}));
expect(eq(fs_suffix(".a."), std::string{"."}));
expect(eq(fs_suffix(".a.b"), std::string{".b"}));
expect(eq(fs_suffix("./b.c"), std::string{".c"}));
expect(eq(fs_suffix("../.b.c"), std::string{".c"}));
expect(eq(fs_suffix("日本語.語"), std::string{".語"}));
expect(eq(fs_suffix("some space.txt"), std::string{".txt"}));

if(fs_is_windows()){
expect(eq(fs_suffix(R"(C:\a\ball.text)"), std::string{".text"}));

if(fs_win32_long_paths_enabled()) {

expect(eq(fs_suffix(R"(\\?\)"), std::string{""}));
expect(eq(fs_suffix(R"(\\.\)"), std::string{""}));

expect(eq(fs_suffix(R"(\\?\C:\)"), std::string{""}));
expect(eq(fs_suffix(R"(\\.\C:\)"), std::string{""}));
expect(eq(fs_suffix(R"(\\?\UNC\server\share)"), std::string{""}));
expect(eq(fs_suffix(R"(\\?\UNC\server\share\日本語.txt)"), std::string{".txt"}));

expect(eq(fs_suffix(R"(\\server\share\some space here.txt)"), std::string{".txt"}));
expect(eq(fs_suffix(R"(\\?\C:\some space here.txt)"), std::string{".txt"}));

expect(fs_win32_is_ext_path(R"(\\.\C:\)"));
expect(fs_win32_is_ext_path(R"(\\?\C:\)"));
}
}

};

}
